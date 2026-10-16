#!/bin/bash
#SBATCH --job-name=mycocosm_download
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=01_genomes/download_%j.log
#SBATCH --error=01_genomes/download_%j.err
# download_mycocosm_genomes.sbatch
# Description: Slurm job that downloads the genomes in the portal list with
#              download_mycocosm_genomes.sh and then decompresses the proteomes
#              with extract_proteomes.sh.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-14
# Version: 1.0
#
# Usage (submit from the project folder): sbatch scripts/download_mycocosm_genomes.sbatch [portal_list]
#   portal_list  Passed on to download_mycocosm_genomes.sh (default: 00_inputs/mycocosm_datasets.tsv).
#
# The Slurm account comes from config.sh: run 'source config.sh' before sbatch. The log is written
# to 01_genomes/download_<jobid>.log. Re-submitting skips files that were already downloaded.

# Slurm starts the job in the folder it was submitted from, which must be the project folder
if [[ ! -d scripts || ! -d 00_inputs ]]; then
    echo "Submit this job from the project folder" >&2
    exit 1
fi

echo "Job $SLURM_JOB_ID on $(hostname), started $(date)"

scripts/download_mycocosm_genomes.sh "$@"
DOWNLOAD_STATUS=$?

# Decompress whatever was downloaded, even if some downloads failed
scripts/extract_proteomes.sh
EXTRACT_STATUS=$?

echo "Finished $(date)"
(( DOWNLOAD_STATUS == 0 && EXTRACT_STATUS == 0 ))
