#!/bin/bash
#SBATCH --job-name=af3
#SBATCH --output=08_structures/logs/af3_%j.out
#SBATCH --error=08_structures/logs/af3_%j.err
#SBATCH --open-mode=append
#SBATCH --partition=gputest
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=217086
#SBATCH --gres=gpu:gh200:1
# run_af3.sbatch
# Description: Slurm job that predicts structures with AlphaFold3 (CSC container, GH200 GPU node)
#              for every input JSON in a folder that has no finished model yet, in one AlphaFold3
#              run (database search and prediction). Based on CSC's template
#              /appl/soft/manual/bio/aarch64/alphafold/3.04/af3_batch_template.sh.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-22
# Version: 1.0
#
# Usage (from the project folder): sbatch [slurm options] scripts/run_af3.sbatch [input_dir] [output_dir]
#   input_dir   AlphaFold3 JSONs from make_af3_inputs.py (default: 08_structures/af3_inputs).
#   output_dir  One subfolder per protein (default: 08_structures/af3_outputs).
# The defaults above (gputest, 15 min) are for tests; for the full set, override them, e.g.
#   sbatch --partition=<gpu partition> --time=<hours> scripts/run_af3.sbatch
#
# A protein is finished when <output_dir>/<name>/<name>_model.cif exists; finished proteins are
# skipped, and the folders of unfinished ones are removed first (AlphaFold3 would otherwise write
# a second, timestamped folder). So a job that runs out of time can simply be resubmitted.
#
# The Slurm account, AF3_MODEL_DIR (folder with af3.bin.zst), AF3_SIF and AF3_DB_DIR come from
# config.sh: run 'source config.sh' before sbatch.

# Slurm starts the job in the folder it was submitted from, which must be the project folder
if [[ ! -d scripts || ! -d 00_inputs ]]; then
    echo "Submit this job from the project folder" >&2
    exit 1
fi

# shellcheck source=config.example.sh
source config.sh || exit 1
for VAR in AF3_MODEL_DIR AF3_SIF AF3_DB_DIR; do
    if [[ -z "${!VAR}" || ! -e "${!VAR}" ]]; then
        echo "$VAR in config.sh is not set or does not exist" >&2
        exit 1
    fi
done
if ! compgen -G "$AF3_MODEL_DIR/af3.bin*" > /dev/null; then
    echo "No af3.bin or af3.bin.zst in AF3_MODEL_DIR ($AF3_MODEL_DIR)" >&2
    exit 1
fi

INPUT_DIR=${1:-08_structures/af3_inputs}
OUTPUT_DIR=${2:-08_structures/af3_outputs}
QUEUE_DIR="08_structures/af3_queue_${SLURM_JOB_ID}"

if ! compgen -G "$INPUT_DIR/*.json" > /dev/null; then
    echo "No JSONs in '$INPUT_DIR', run scripts/make_af3_inputs.py first" >&2
    exit 1
fi

# Queue: links to the inputs without a finished model; unfinished output folders are removed
mkdir -p "$OUTPUT_DIR" "$QUEUE_DIR"
N_INPUTS=0
N_QUEUED=0
for JSON in "$INPUT_DIR"/*.json; do
    NAME=$(basename "$JSON" .json)
    N_INPUTS=$(( N_INPUTS + 1 ))
    [[ -s "$OUTPUT_DIR/$NAME/${NAME}_model.cif" ]] && continue
    rm -rf "${OUTPUT_DIR:?}/$NAME"
    ln -s "$(realpath "$JSON")" "$QUEUE_DIR/$NAME.json"
    N_QUEUED=$(( N_QUEUED + 1 ))
done

echo "Job $SLURM_JOB_ID on $(hostname), started $(date)"
echo "$N_QUEUED of $N_INPUTS proteins to predict from $INPUT_DIR into $OUTPUT_DIR"
if (( N_QUEUED > 0 )); then
    export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
    apptainer exec --nv \
        --env "PATH=/hmmer/bin:/opt/alphafold3/.venv/bin:$PATH" \
        --bind "$(csc-common-bind)" \
        "$AF3_SIF" \
        /opt/alphafold3/.venv/bin/python /opt/alphafold3/run_alphafold.py \
            --input_dir="$QUEUE_DIR" \
            --model_dir="$AF3_MODEL_DIR" \
            --db_dir="$AF3_DB_DIR" \
            --output_dir="$OUTPUT_DIR"
    echo "AlphaFold3 finished $(date) with exit code $?"
fi
rm -rf "$QUEUE_DIR"

N_DONE=0
for JSON in "$INPUT_DIR"/*.json; do
    NAME=$(basename "$JSON" .json)
    [[ -s "$OUTPUT_DIR/$NAME/${NAME}_model.cif" ]] && N_DONE=$(( N_DONE + 1 ))
done
echo "Finished $(date): $N_DONE of $N_INPUTS proteins have a model"
if (( N_DONE < N_INPUTS )); then
    echo "$(( N_INPUTS - N_DONE )) proteins have no model yet; resubmit to continue" >&2
    exit 1
fi
