#!/bin/bash
#SBATCH --job-name=family_search
#SBATCH --output=09_families/logs/search_%j.out
#SBATCH --error=09_families/logs/search_%j.err
#SBATCH --open-mode=append
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
# search_family.sbatch
# Description: Slurm job with the two searches behind a gene family (families/<family>/):
#              1. DIAMOND blastp (ultra-sensitive) of the family's reference proteins
#                 (references.tsv, use = TRUE) against all cohort proteins, to find their copies
#                 (anchors) in the cohort;
#              2. hmmsearch of each of the family's Pfam HMMs (HMMS in family.conf) against all
#                 cohort proteins and against the references, at Pfam's gathering thresholds.
#              summarize_family.py turns the results into anchors, orphans, members and orthogroups.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-23
# Version: 1.0
#
# Usage (from the project folder): sbatch scripts/search_family.sbatch <family> [proteomes_dir]
#   family         Folder name in families/ (e.g. sugar_transporters).
#   proteomes_dir  Clean cohort proteomes (default: 05_clean_proteomes/cohort/proteomes).
#
# Output, inside <OUT_DIR>/search (OUT_DIR from family.conf):
#   reference_proteins.faa        the references searched
#   diamond_references.tsv        DIAMOND hits: qseqid sseqid pident length qlen slen qstart qend
#                                 sstart send evalue bitscore (up to 50 targets per reference)
#   <HMM>.hmm                     each HMM, downloaded from InterPro if missing (if the compute node
#                                 has no internet, run the curl command printed in the log on the
#                                 login node)
#   <hmm>_cohort.domtbl, <hmm>_references.domtbl   hmmsearch per-domain tables
# Finished searches are skipped; if the references change, their searches are repeated.
#
# The Slurm account comes from config.sh: run 'source config.sh' before sbatch.

# Slurm starts the job in the folder it was submitted from, which must be the project folder
if [[ ! -d scripts || ! -d 00_inputs ]]; then
    echo "Submit this job from the project folder" >&2
    exit 1
fi

FAMILY=$1
PROTEOMES_DIR=${2:-05_clean_proteomes/cohort/proteomes}
CONF="families/$FAMILY/family.conf"
REFERENCES="families/$FAMILY/references.tsv"
if [[ -z "$FAMILY" || ! -s "$CONF" || ! -s "$REFERENCES" ]]; then
    echo "Usage: sbatch scripts/search_family.sbatch <family>; families/<family>/ needs family.conf and references.tsv" >&2
    exit 1
fi
# shellcheck source=/dev/null
source "$CONF"
OUT="$OUT_DIR/search"
THREADS=${SLURM_CPUS_PER_TASK:-4}

if ! compgen -G "$PROTEOMES_DIR/*.faa" > /dev/null; then
    echo "No proteomes in '$PROTEOMES_DIR', run scripts/clean_proteomes.py first" >&2
    exit 1
fi

module load bio-apps/v202603
module load diamond/2.2.5
module load hmmer/3.4
for TOOL in diamond hmmsearch; do
    if ! command -v "$TOOL" > /dev/null; then
        echo "$TOOL not found, check that the modules above loaded" >&2
        exit 1
    fi
done

set -o pipefail
mkdir -p "$OUT"
echo "Job $SLURM_JOB_ID on $(hostname), started $(date): $NAME ($HMMS); $(diamond version | head -1)," \
    "$(hmmsearch -h | sed -n 2p | sed 's/^# //')"

# References in use; when they changed since the last run, their searches are repeated
awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
    $col["use"] == "TRUE" { print ">" $col["reference_id"]; print $col["sequence"] }' \
    "$REFERENCES" | tr -d '\r' > "$OUT/reference_proteins.faa.new"
if ! cmp -s "$OUT/reference_proteins.faa.new" "$OUT/reference_proteins.faa"; then
    rm -f "$OUT/diamond_references.tsv" "$OUT"/*_references.domtbl
fi
mv "$OUT/reference_proteins.faa.new" "$OUT/reference_proteins.faa"
echo "$(grep -c '^>' "$OUT/reference_proteins.faa") reference proteins"

# All cohort proteins in one FASTA, removed at the end
ALL="$OUT/cohort_proteins.faa"
cat "$PROTEOMES_DIR"/*.faa > "$ALL" || exit 1
echo "$(grep -c '^>' "$ALL") cohort proteins from $(find "$PROTEOMES_DIR" -maxdepth 1 -name '*.faa' | wc -l) proteomes"

# 1. DIAMOND search of the references
if [[ -s "$OUT/diamond_references.tsv" ]]; then
    echo "DIAMOND results already exist"
else
    diamond makedb --in "$ALL" -d "$OUT/cohort_proteins" --threads "$THREADS" --quiet || exit 1
    if ! diamond blastp -q "$OUT/reference_proteins.faa" -d "$OUT/cohort_proteins" --ultra-sensitive \
        --max-target-seqs 50 --threads "$THREADS" --quiet -o "$OUT/diamond_references.tsv.part" \
        --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore; then
        echo "DIAMOND failed" >&2
        exit 1
    fi
    mv "$OUT/diamond_references.tsv.part" "$OUT/diamond_references.tsv"
    rm -f "$OUT/cohort_proteins.dmnd"
    echo "DIAMOND: $(cut -f1 "$OUT/diamond_references.tsv" | sort -u | wc -l) references with hits"
fi

# 2. HMM searches
for HMM in $HMMS; do
    hmm=$(echo "$HMM" | tr '[:upper:]' '[:lower:]')
    URL="https://www.ebi.ac.uk/interpro/api/entry/pfam/$HMM?annotation=hmm"
    if [[ ! -s "$OUT/$HMM.hmm" ]]; then
        if ! curl -sfL "$URL" | gunzip > "$OUT/$HMM.hmm.part" || ! grep -q "^HMMER3" "$OUT/$HMM.hmm.part"; then
            rm -f "$OUT/$HMM.hmm.part"
            echo "Could not download $HMM; on the login node run:" >&2
            echo "  curl -sfL '$URL' | gunzip > $OUT/$HMM.hmm" >&2
            exit 1
        fi
        mv "$OUT/$HMM.hmm.part" "$OUT/$HMM.hmm"
    fi
    echo "HMM $(grep -m1 '^ACC' "$OUT/$HMM.hmm" | awk '{print $2}'), length $(grep -m1 '^LENG' "$OUT/$HMM.hmm" | awk '{print $2}')"
    for TARGET in cohort references; do
        TABLE="$OUT/${hmm}_$TARGET.domtbl"
        [[ -s "$TABLE" ]] && { echo "$TABLE already exists"; continue; }
        INPUT=$([[ $TARGET == cohort ]] && echo "$ALL" || echo "$OUT/reference_proteins.faa")
        if ! hmmsearch --cut_ga --cpu "$THREADS" --noali -o /dev/null --domtblout "$TABLE.part" \
            "$OUT/$HMM.hmm" "$INPUT"; then
            echo "hmmsearch of $HMM failed on $TARGET" >&2
            rm -f "$TABLE.part"
            exit 1
        fi
        # hmmsearch records the working directory in a comment line; it is a local path, so it is dropped
        grep -v '^# Current dir:' "$TABLE.part" > "$TABLE" && rm "$TABLE.part"
        echo "$HMM on $TARGET: $(grep -v '^#' "$TABLE" | awk '{print $1}' | sort -u | wc -l) proteins with hits"
    done
done

rm -f "$ALL"
echo "Done $(date). Next: python3 scripts/summarize_family.py $FAMILY"
