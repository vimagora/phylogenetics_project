#!/bin/bash
#SBATCH --job-name=orthofinder
#SBATCH --output=07_orthofinder/logs/orthofinder_%j.out
#SBATCH --error=07_orthofinder/logs/orthofinder_%j.err
#SBATCH --open-mode=append
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=256G
# run_orthofinder.sbatch
# Description: Slurm job that runs OrthoFinder 3 on the clean cohort proteomes with our rooted
#              species tree (-s), the default DIAMOND search and MSA gene trees (FAMSA + FastTree).
#              Sequence IDs are kept as they are (-X), so genes stay '<portal>-<proteinId>'.
#              When OrthoFinder finishes, the folders with one file per orthogroup or species
#              pair are packed into one .tar.gz each, to keep the scratch file count low.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-20
# Version: 1.0
#
# Usage (from the project folder): sbatch scripts/run_orthofinder.sbatch [proteomes_dir] [species_tree]
#   proteomes_dir  One FASTA per species; the file name is the species name
#                  (default: 05_clean_proteomes/cohort/proteomes).
#   species_tree   Rooted topology from prepare_species_tree.py, tips = file names without .faa
#                  (default: 07_orthofinder/species_tree_rooted.nwk).
#
# Results go to 07_orthofinder/run/Results_cohort/. OrthoFinder refuses an existing output folder,
# so the job does not start if 07_orthofinder/run exists; OrthoFinder cannot resume a run. If a run stops after the DIAMOND search, its
# search results can be reused with 'orthofinder -b <Results>/WorkingDirectory' (see README).
#
# Packed after the run (extract single files with 'tar -xzf <folder>.tar.gz <folder>/<file>'):
#   WorkingDirectory, Orthogroup_Sequences, MultipleSequenceAlignments, Gene_Trees,
#   Resolved_Gene_Trees, Orthologues
#
# The Slurm account and ORTHOFINDER_DIR come from config.sh: run 'source config.sh' before sbatch.

# Slurm starts the job in the folder it was submitted from, which must be the project folder
if [[ ! -d scripts || ! -d 00_inputs ]]; then
    echo "Submit this job from the project folder" >&2
    exit 1
fi

# Local settings: ORTHOFINDER_DIR, the Tykky container with OrthoFinder 3
# shellcheck source=config.example.sh
source config.sh || exit 1
if [[ -z "$ORTHOFINDER_DIR" || ! -x "$ORTHOFINDER_DIR/bin/orthofinder" ]]; then
    echo "ORTHOFINDER_DIR in config.sh does not contain bin/orthofinder" >&2
    exit 1
fi
export PATH="$ORTHOFINDER_DIR/bin:$PATH"

# The bioconda launcher imports its code from a 'src' folder next to itself, but the Tykky wrapper
# runs it under the wrapper's path, where there is none. Point Python to the folder inside the
# container instead; Tykky's mount point (with a random suffix) is in common.sh.
TYKKY_ROOT=$(sed -n 's/^INSTALLATION_PATH=//p' "$ORTHOFINDER_DIR/common.sh" 2> /dev/null)
if [[ -z "$TYKKY_ROOT" ]]; then
    echo "No INSTALLATION_PATH in $ORTHOFINDER_DIR/common.sh; is ORTHOFINDER_DIR a Tykky container?" >&2
    exit 1
fi
export PYTHONPATH="$TYKKY_ROOT/miniforge/envs/env1/bin/src"
if ! orthofinder -v > /dev/null 2>&1; then
    echo "orthofinder does not start (PYTHONPATH=$PYTHONPATH); test it with 'orthofinder -v'" >&2
    exit 1
fi

PROTEOMES_DIR=${1:-05_clean_proteomes/cohort/proteomes}
SPECIES_TREE=${2:-07_orthofinder/species_tree_rooted.nwk}
RUN_DIR=07_orthofinder/run
RUN_NAME=cohort
RESULTS="$RUN_DIR/Results_$RUN_NAME"

if ! compgen -G "$PROTEOMES_DIR/*.faa" > /dev/null; then
    echo "No proteomes in '$PROTEOMES_DIR', run scripts/clean_proteomes.py first" >&2
    exit 1
fi
if [[ ! -s "$SPECIES_TREE" ]]; then
    echo "Species tree '$SPECIES_TREE' not found, run scripts/prepare_species_tree.py first" >&2
    exit 1
fi
if [[ -e "$RUN_DIR" ]]; then
    echo "'$RUN_DIR' already exists; move or delete it to start a new run" >&2
    exit 1
fi

# OrthoFinder opens many files at once (one per species pair during the search)
ulimit -n "$(ulimit -Hn)"

echo "Job $SLURM_JOB_ID on $(hostname), started $(date)"
echo "OrthoFinder $(orthofinder -v 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1) on" \
    "$(find "$PROTEOMES_DIR" -maxdepth 1 -name '*.faa' | wc -l) proteomes in $PROTEOMES_DIR," \
    "species tree $SPECIES_TREE"

orthofinder -f "$PROTEOMES_DIR" -s "$SPECIES_TREE" -X -t "$SLURM_CPUS_PER_TASK" -a 16 \
    -o "$RUN_DIR" -n "$RUN_NAME"
OF_STATUS=$?
echo "OrthoFinder finished $(date) with exit code $OF_STATUS"
# OrthoFinder can exit with 0 after an error (e.g. disk quota exceeded), and Orthogroups.tsv is
# written early, so success is judged by the last line OrthoFinder writes to its Log.txt
if (( OF_STATUS != 0 )) || ! grep -q "OrthoFinder run completed" "$RESULTS/Log.txt" 2> /dev/null; then
    echo "OrthoFinder did not complete, see the log above and $RESULTS/Log.txt;" \
        "the results are incomplete and nothing was packed" >&2
    exit 1
fi

# Pack the folders with one file per orthogroup or species pair; each is only removed once
# its archive is complete
for FOLDER in WorkingDirectory Orthogroup_Sequences MultipleSequenceAlignments Gene_Trees \
    Resolved_Gene_Trees Orthologues; do
    [[ -d "$RESULTS/$FOLDER" ]] || continue
    N_FILES=$(find "$RESULTS/$FOLDER" -type f | wc -l)
    if tar -czf "$RESULTS/$FOLDER.tar.gz.part" -C "$RESULTS" "$FOLDER" \
        && mv "$RESULTS/$FOLDER.tar.gz.part" "$RESULTS/$FOLDER.tar.gz"; then
        rm -r "${RESULTS:?}/$FOLDER"
        echo "Packed $FOLDER ($N_FILES files) into $RESULTS/$FOLDER.tar.gz"
    else
        echo "Could not pack $FOLDER, left as it is" >&2
        rm -f "$RESULTS/$FOLDER.tar.gz.part"
    fi
done

echo "Done $(date). Orthogroups: $RESULTS/Orthogroups/Orthogroups.tsv"
