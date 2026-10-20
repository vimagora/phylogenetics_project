# shellcheck shell=bash disable=SC2034
# config.sh: local settings that stay out of git.
# Copy this file to config.sh, fill in your values, and run 'source config.sh' before sbatch.

# Slurm account (project) for sbatch; sbatch reads it from the environment
export SBATCH_ACCOUNT=project_XXXXXXX

# Folder with the BUSCO datasets, i.e. the folder that contains lineages/fungi_odb12
BUSCO_DOWNLOADS=/path/to/busco_downloads/

# Tykky container with OrthoFinder 3 (envs/of3_env.yaml), i.e. the folder that contains bin/orthofinder
ORTHOFINDER_DIR=/path/to/orthofinder/
