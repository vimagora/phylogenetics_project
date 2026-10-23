#!/bin/bash
#SBATCH --job-name=family_align
#SBATCH --output=09_families/logs/align_%j.out
#SBATCH --error=09_families/logs/align_%j.err
#SBATCH --open-mode=append
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
# align_family.sbatch
# Description: Slurm job that aligns a gene family's tree sequences (build_family_fasta.py) with
#              MAFFT --auto, which picks a strategy for the number of sequences (FFT-NS-2 for
#              thousands), and trims the alignment with ClipKIT in gappy mode (columns with more than
#              50% gaps removed), as for the species tree. Columns of rare extra domains, e.g. from
#              fusions, are mostly gaps and are trimmed away.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-23
# Version: 1.0
#
# Usage (from the project folder): sbatch scripts/align_family.sbatch <family>
#
# Output, inside <OUT_DIR>/tree (OUT_DIR from families/<family>/family.conf):
#   aligned.faa   MAFFT alignment (the strategy MAFFT chose is in the job's .err log)
#   trimmed.faa   ClipKIT-trimmed alignment, the input of family_tree.sbatch
# An existing alignment is reused; delete it to realign.
#
# The Slurm account comes from config.sh: run 'source config.sh' before sbatch.

# Slurm starts the job in the folder it was submitted from, which must be the project folder
if [[ ! -d scripts || ! -d 00_inputs ]]; then
    echo "Submit this job from the project folder" >&2
    exit 1
fi

FAMILY=$1
CONF="families/$FAMILY/family.conf"
if [[ -z "$FAMILY" || ! -s "$CONF" ]]; then
    echo "Usage: sbatch scripts/align_family.sbatch <family>" >&2
    exit 1
fi
# shellcheck source=/dev/null
source "$CONF"
DIR="$OUT_DIR/tree"
TRIM_GAPS=0.5
if [[ ! -s "$DIR/sequences.faa" ]]; then
    echo "'$DIR/sequences.faa' not found, run scripts/build_family_fasta.py $FAMILY first" >&2
    exit 1
fi

module load bio-apps/v202603
module load mafft/7.525
module load clipkit/2.12.1
for TOOL in mafft clipkit; do
    if ! command -v "$TOOL" > /dev/null; then
        echo "$TOOL not found, check that the modules above loaded" >&2
        exit 1
    fi
done

echo "Job $SLURM_JOB_ID on $(hostname), started $(date): $(grep -c '^>' "$DIR/sequences.faa") sequences of $NAME"
if [[ -s "$DIR/aligned.faa" ]]; then
    echo "Alignment already exists: $DIR/aligned.faa"
else
    if ! mafft --auto --anysymbol --thread "$SLURM_CPUS_PER_TASK" "$DIR/sequences.faa" > "$DIR/aligned.faa.part"; then
        echo "MAFFT failed" >&2
        rm -f "$DIR/aligned.faa.part"
        exit 1
    fi
    mv "$DIR/aligned.faa.part" "$DIR/aligned.faa"
fi
if ! clipkit "$DIR/aligned.faa" -m gappy -g "$TRIM_GAPS" -of fasta -o "$DIR/trimmed.faa.part" > /dev/null; then
    echo "ClipKIT failed" >&2
    exit 1
fi
mv "$DIR/trimmed.faa.part" "$DIR/trimmed.faa"
length() { awk '/^>/ { if (seen++) exit; next } { n += length($0) } END { print n }' "$1"; }
echo "Done $(date): $(length "$DIR/aligned.faa") aligned columns, $(length "$DIR/trimmed.faa") after trimming." \
    "Next: sbatch --time=04:00:00 scripts/family_tree.sbatch $FAMILY fasttree"
