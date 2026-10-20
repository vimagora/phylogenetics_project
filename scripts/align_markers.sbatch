#!/bin/bash
#SBATCH --job-name=align_markers
#SBATCH --output=06_species_tree/logs/align_%j.out
#SBATCH --error=06_species_tree/logs/align_%j.err
#SBATCH --open-mode=append
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
# align_markers.sbatch
# Description: Slurm job that aligns every marker FASTA with MAFFT L-INS-i and trims the alignment
#              with ClipKIT (gappy mode: columns with more than TRIM_GAPS gaps are removed), running
#              one marker per CPU in parallel. ClipKIT's default smart-gap mode kept most columns
#              of these kingdom-wide alignments (median 46% gaps per trimmed marker).
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-20
# Version: 1.0
#
# Usage (from the project folder): sbatch scripts/align_markers.sbatch [markers_dir]
#   markers_dir  Folder with raw/<busco_id>.faa from extract_marker_sequences.py
#                (default: 06_species_tree/markers).
#
# Output, inside markers_dir:
#   aligned/<busco_id>.faa  MAFFT alignment
#   trimmed/<busco_id>.faa  ClipKIT-trimmed alignment
# Markers already aligned or trimmed are skipped, so a failed or requeued job can be resubmitted.
# To re-trim with other settings, delete the trimmed folder and resubmit: alignments are reused.
#
# The Slurm account comes from config.sh: run 'source config.sh' before sbatch.

# Slurm starts the job in the folder it was submitted from, which must be the project folder
if [[ ! -d scripts || ! -d 00_inputs ]]; then
    echo "Submit this job from the project folder" >&2
    exit 1
fi

MARKERS_DIR=${1:-06_species_tree/markers}
RAW_DIR="$MARKERS_DIR/raw"
ALN_DIR="$MARKERS_DIR/aligned"
TRIM_DIR="$MARKERS_DIR/trimmed"

# ClipKIT gappy mode: columns with a gap fraction above TRIM_GAPS are removed
TRIM_GAPS=0.5

if ! compgen -G "$RAW_DIR/*.faa" > /dev/null; then
    echo "No marker FASTAs in '$RAW_DIR', run scripts/extract_marker_sequences.py first" >&2
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

echo "Job $SLURM_JOB_ID on $(hostname), started $(date)"
echo "Aligning and trimming $(find "$RAW_DIR" -maxdepth 1 -name '*.faa' | wc -l) markers in $RAW_DIR" \
    "(ClipKIT gappy mode, gap threshold $TRIM_GAPS)"
mkdir -p "$ALN_DIR" "$TRIM_DIR"

# Aligns and trims one marker. Outputs go to .part files that are only renamed when the tool
# succeeds, so an interrupted job never leaves a truncated alignment that would be skipped later.
align_one() {
    local id=$1
    local raw="$RAW_DIR/$id.faa" aln="$ALN_DIR/$id.faa" trim="$TRIM_DIR/$id.faa"
    if [[ ! -s "$aln" ]]; then
        if ! mafft --localpair --maxiterate 1000 --thread 1 --quiet "$raw" > "$aln.part"; then
            echo "  FAILED mafft: $id" >&2
            rm -f "$aln.part"
            return 1
        fi
        mv "$aln.part" "$aln"
    fi
    if [[ ! -s "$trim" ]]; then
        if ! clipkit "$aln" -m gappy -g "$TRIM_GAPS" -of fasta -o "$trim.part" > /dev/null; then
            echo "  FAILED clipkit: $id" >&2
            rm -f "$trim.part"
            return 1
        fi
        mv "$trim.part" "$trim"
    fi
}
export -f align_one
export RAW_DIR ALN_DIR TRIM_DIR TRIM_GAPS

find "$RAW_DIR" -maxdepth 1 -name '*.faa' -printf '%f\n' | sed 's/\.faa$//' | sort \
    | xargs -P "$SLURM_CPUS_PER_TASK" -I {} bash -c 'align_one "$1"' _ {}

N_RAW=$(find "$RAW_DIR" -maxdepth 1 -name '*.faa' | wc -l)
N_TRIM=$(find "$TRIM_DIR" -maxdepth 1 -name '*.faa' | wc -l)
echo "Finished $(date): $N_TRIM of $N_RAW markers aligned and trimmed"
if (( N_TRIM < N_RAW )); then
    echo "$(( N_RAW - N_TRIM )) markers failed, see the errors above; resubmit to retry them" >&2
    exit 1
fi
