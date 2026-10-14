#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --account=project_2015320
#SBATCH --output=/scratch/project_2015320/phylogenetics_project/02_BUSCO/logs/busco_%j.out
#SBATCH --error=/scratch/project_2015320/phylogenetics_project/02_BUSCO/logs/busco_%j.err
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
# run_busco.sbatch
# Description: Slurm job that runs BUSCO in batch mode (protein mode, fungi_odb12) on all the
#              decompressed proteomes in a folder, one after the other.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-14
# Version: 1.0
#
# Usage (from the project folder): sbatch scripts/run_busco.sbatch [proteomes_dir]
#   proteomes_dir  Folder with the decompressed proteomes (default: 02_BUSCO/proteomes).
#
# Results go to 02_BUSCO/results/<proteomes_dir name>/, with one subfolder per proteome and a
# batch_summary.txt. BUSCO refuses to start if that results folder already exists.

cd /scratch/project_2015320/phylogenetics_project/ || exit 1

PROTEOMES_DIR=${1:-02_BUSCO/proteomes}
RESULTS_DIR=02_BUSCO/results
BUSCO_DOWNLOADS=/projappl/project_2015320/busco_downloads/

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
