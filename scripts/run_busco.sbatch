#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --output=03_BUSCO/logs/busco_%j.out
#SBATCH --error=03_BUSCO/logs/busco_%j.err
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
# run_busco.sbatch
# Description: Slurm job that runs BUSCO in batch mode (protein mode, fungi_odb12) on all the
#              decompressed proteomes in a folder, one after the other.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-14
# Version: 1.0
#
# Usage (from the project folder): sbatch scripts/run_busco.sbatch [proteomes_dir]
#   proteomes_dir  Folder with the decompressed proteomes (default: 03_BUSCO/proteomes).
#
# Results go to 03_BUSCO/results/<proteomes_dir name>/, with one subfolder per proteome and a
# batch_summary.txt. BUSCO refuses to start if that results folder already exists.
#
# The Slurm account and the BUSCO datasets folder come from config.sh: run 'source config.sh'
# before sbatch.

# Slurm starts the job in the folder it was submitted from, which must be the project folder
if [[ ! -d scripts || ! -d 00_inputs ]]; then
    echo "Submit this job from the project folder" >&2
    exit 1
fi

# Local settings: BUSCO_DOWNLOADS, the folder with the BUSCO datasets (lineages/fungi_odb12)
# shellcheck source=config.example.sh
source config.sh || exit 1
if [[ -z "$BUSCO_DOWNLOADS" ]]; then
    echo "BUSCO_DOWNLOADS is not set in config.sh" >&2
    exit 1
fi

PROTEOMES_DIR=${1:-03_BUSCO/proteomes}
RESULTS_DIR=03_BUSCO/results

if [[ ! -d "$PROTEOMES_DIR" ]]; then
    echo "Proteomes folder '$PROTEOMES_DIR' not found, run scripts/extract_proteomes.sh first" >&2
    exit 1
fi

# BUSCO runs on every file in the folder, so leftovers of an interrupted extraction would be included
if compgen -G "$PROTEOMES_DIR/*.part" > /dev/null; then
    echo "Unfinished .part files in '$PROTEOMES_DIR', re-run scripts/extract_proteomes.sh first" >&2
    exit 1
fi

echo "Job $SLURM_JOB_ID on $(hostname), started $(date)"
echo "Running BUSCO on $(find "$PROTEOMES_DIR" -maxdepth 1 -type f | wc -l) files in $PROTEOMES_DIR"

module load bio-apps/v202603
module load busco/6.1.0

mkdir -p "$RESULTS_DIR"
busco -i "$PROTEOMES_DIR" -m proteins -l fungi_odb12 --download_path "$BUSCO_DOWNLOADS" \
    --out_path "$RESULTS_DIR" -o "$(basename "$PROTEOMES_DIR")" -c "$SLURM_CPUS_PER_TASK"
BUSCO_STATUS=$?

echo "Finished $(date) with exit code $BUSCO_STATUS"
exit $BUSCO_STATUS
