#!/bin/bash
#SBATCH --job-name=st_search
#SBATCH --output=09_sugar_transporters/logs/st_search_%j.out
#SBATCH --error=09_sugar_transporters/logs/st_search_%j.err
#SBATCH --open-mode=append
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
# search_sugar_transporters.sbatch
# Description: Slurm job with the two searches behind the cohort sugar transporters:
#              1. DIAMOND blastp (ultra-sensitive) of the reference transporters (characterized and
#                 D. squalens, from 08_structures/structure_ids.tsv) against all cohort proteins,
#                 to find their copies (anchors) in the cohort;
#              2. hmmsearch of Pfam PF00083 (Sugar_tr) against all cohort proteins and against the
#                 references, with Pfam's gathering thresholds (--cut_ga).
#              summarize_sugar_transporters.py turns the results into anchors, sugar porters and
#              orthogroups.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-23
# Version: 1.0
#
# Usage (from the project folder): sbatch scripts/search_sugar_transporters.sbatch [proteomes_dir]
#   proteomes_dir  Clean cohort proteomes (default: 05_clean_proteomes/cohort/proteomes).
#
# Output, inside 09_sugar_transporters/search:
#   reference_proteins.faa      the references searched (use_in_structure_analysis = TRUE)
#   diamond_references.tsv      DIAMOND hits: qseqid sseqid pident length qlen slen qstart qend
#                               sstart send evalue bitscore (up to 50 targets per reference)
#   PF00083.hmm                 the Pfam HMM, downloaded from InterPro if missing (if the compute
#                               node has no internet, download it on the login node with the curl
#                               command below)
#   pf00083_cohort.domtbl       hmmsearch per-domain table for the cohort proteins
#   pf00083_references.domtbl   the same for the references
# Finished searches are skipped, so the job can be resubmitted.
#
# The Slurm account comes from config.sh: run 'source config.sh' before sbatch.

# Slurm starts the job in the folder it was submitted from, which must be the project folder
if [[ ! -d scripts || ! -d 00_inputs ]]; then
    echo "Submit this job from the project folder" >&2
    exit 1
fi

PROTEOMES_DIR=${1:-05_clean_proteomes/cohort/proteomes}
IDS_TABLE=08_structures/structure_ids.tsv
OUT_DIR=09_sugar_transporters/search
HMM_URL="https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00083?annotation=hmm"
THREADS=${SLURM_CPUS_PER_TASK:-4}

if ! compgen -G "$PROTEOMES_DIR/*.faa" > /dev/null; then
    echo "No proteomes in '$PROTEOMES_DIR', run scripts/clean_proteomes.py first" >&2
    exit 1
fi
if [[ ! -s "$IDS_TABLE" ]]; then
    echo "'$IDS_TABLE' not found, run scripts/make_af3_inputs.py first" >&2
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
mkdir -p "$OUT_DIR"
echo "Job $SLURM_JOB_ID on $(hostname), started $(date); $(diamond version | head -1), $(hmmsearch -h | sed -n 2p)"

# References: the proteins used in the structure analysis, with their full sequence
awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
    $col["use_in_structure_analysis"] == "TRUE" { print ">" $col["structure_id"]; print $col["sequence"] }' \
    "$IDS_TABLE" > "$OUT_DIR/reference_proteins.faa"
echo "$(grep -c '^>' "$OUT_DIR/reference_proteins.faa") reference proteins"

# All cohort proteins in one FASTA, removed at the end
ALL="$OUT_DIR/cohort_proteins.faa"
cat "$PROTEOMES_DIR"/*.faa > "$ALL" || exit 1
echo "$(grep -c '^>' "$ALL") cohort proteins from $(find "$PROTEOMES_DIR" -maxdepth 1 -name '*.faa' | wc -l) proteomes"

# 1. DIAMOND search of the references
if [[ -s "$OUT_DIR/diamond_references.tsv" ]]; then
    echo "DIAMOND results already exist"
else
    diamond makedb --in "$ALL" -d "$OUT_DIR/cohort_proteins" --threads "$THREADS" --quiet || exit 1
    if ! diamond blastp -q "$OUT_DIR/reference_proteins.faa" -d "$OUT_DIR/cohort_proteins" --ultra-sensitive \
        --max-target-seqs 50 --threads "$THREADS" --quiet -o "$OUT_DIR/diamond_references.tsv.part" \
        --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore; then
        echo "DIAMOND failed" >&2
        exit 1
    fi
    mv "$OUT_DIR/diamond_references.tsv.part" "$OUT_DIR/diamond_references.tsv"
    rm -f "$OUT_DIR/cohort_proteins.dmnd"
    echo "DIAMOND: $(cut -f1 "$OUT_DIR/diamond_references.tsv" | sort -u | wc -l) references with hits"
fi

# 2. PF00083 searches
if [[ ! -s "$OUT_DIR/PF00083.hmm" ]]; then
    if ! curl -sfL "$HMM_URL" | gunzip > "$OUT_DIR/PF00083.hmm.part" || ! grep -q "^HMMER3" "$OUT_DIR/PF00083.hmm.part"; then
        rm -f "$OUT_DIR/PF00083.hmm.part"
        echo "Could not download PF00083; on the login node run:" >&2
        echo "  curl -sfL '$HMM_URL' | gunzip > $OUT_DIR/PF00083.hmm" >&2
        exit 1
    fi
    mv "$OUT_DIR/PF00083.hmm.part" "$OUT_DIR/PF00083.hmm"
fi
echo "HMM: $(grep -m1 '^ACC' "$OUT_DIR/PF00083.hmm" | awk '{print $2}'), length $(grep -m1 '^LENG' "$OUT_DIR/PF00083.hmm" | awk '{print $2}')"
for TARGET in cohort references; do
    TABLE="$OUT_DIR/pf00083_$TARGET.domtbl"
    [[ -s "$TABLE" ]] && { echo "$TABLE already exists"; continue; }
    INPUT=$([[ $TARGET == cohort ]] && echo "$ALL" || echo "$OUT_DIR/reference_proteins.faa")
    if ! hmmsearch --cut_ga --cpu "$THREADS" --noali -o /dev/null --domtblout "$TABLE.part" \
        "$OUT_DIR/PF00083.hmm" "$INPUT"; then
        echo "hmmsearch failed on $TARGET" >&2
        rm -f "$TABLE.part"
        exit 1
    fi
    mv "$TABLE.part" "$TABLE"
    echo "PF00083 on $TARGET: $(grep -v '^#' "$TABLE" | awk '{print $1}' | sort -u | wc -l) proteins with hits"
done

rm -f "$ALL"
echo "Done $(date). Next: python3 scripts/summarize_sugar_transporters.py"
