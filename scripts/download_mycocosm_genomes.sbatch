#!/bin/bash
#SBATCH --job-name=mycocosm_download
#SBATCH --account=project_2015320
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/scratch/project_2015320/phylogenetics_project/01_genomes/download_%j.log
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
# Fill in --account and --partition above before submitting. The log is written to
# 01_genomes/download_<jobid>.log. Re-submitting skips files that were already downloaded.

cd /scratch/project_2015320/phylogenetics_project/ || exit 1

echo "Job $SLURM_JOB_ID on $(hostname), started $(date)"

scripts/download_mycocosm_genomes.sh "$@"
DOWNLOAD_STATUS=$?

# Decompress whatever was downloaded, even if some downloads failed
scripts/extract_proteomes.sh
EXTRACT_STATUS=$?

echo "Finished $(date)"
(( DOWNLOAD_STATUS == 0 && EXTRACT_STATUS == 0 ))
