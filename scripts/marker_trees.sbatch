#!/bin/bash
#SBATCH --job-name=marker_trees
#SBATCH --output=06_species_tree/logs/marker_trees_%j.out
#SBATCH --error=06_species_tree/logs/marker_trees_%j.err
#SBATCH --open-mode=append
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=192G
# marker_trees.sbatch
# Description: Slurm job that infers one maximum-likelihood tree per kept marker with IQ-TREE 3
#              (ModelFinder among nuclear models, 1000 ultrafast bootstraps with NNI refinement),
#              running one marker per CPU in parallel. The trees are the input of ASTRAL, and
#              their best-fit models are reused by the concatenated run.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-20
# Version: 1.0
#
# Usage (from the project folder): sbatch scripts/marker_trees.sbatch [markers_list] [trimmed_dir] [output_dir]
#   markers_list  Markers to use (default: 06_species_tree/concatenated/kept_markers.txt).
#   trimmed_dir   Trimmed alignments (default: 06_species_tree/markers/trimmed).
#   output_dir    Where the trees are saved (default: 06_species_tree/marker_trees).
#
# Output, inside output_dir: <busco_id>.treefile (ML tree with bootstrap support), .iqtree (report,
# with the best-fit model) and .log. Other IQ-TREE files are deleted once a marker is finished,
# to save disk space and file count. Finished markers are skipped, so a failed, timed-out or
# requeued job can be resubmitted.
#
# The Slurm account comes from config.sh: run 'source config.sh' before sbatch.

# Slurm starts the job in the folder it was submitted from, which must be the project folder
if [[ ! -d scripts || ! -d 00_inputs ]]; then
    echo "Submit this job from the project folder" >&2
    exit 1
fi

MARKERS_LIST=${1:-06_species_tree/concatenated/kept_markers.txt}
TRIM_DIR=${2:-06_species_tree/markers/trimmed}
TREES_DIR=${3:-06_species_tree/marker_trees}

if [[ ! -s "$MARKERS_LIST" ]]; then
    echo "Markers list '$MARKERS_LIST' not found, run scripts/concatenate_markers.py first" >&2
    exit 1
fi

module load bio-apps/v202603
module load iqtree3/3.1.4
if ! command -v iqtree3 > /dev/null; then
    echo "iqtree3 not found, check that the modules above loaded" >&2
    exit 1
fi

echo "Job $SLURM_JOB_ID on $(hostname), started $(date)"
echo "Inferring trees for $(grep -c . "$MARKERS_LIST") markers from $TRIM_DIR"
mkdir -p "$TREES_DIR"

# Infers one marker's tree. An unfinished run is resumed from its IQ-TREE checkpoint; a marker
# counts as finished once its .iqtree report exists, and its intermediate files are then removed.
tree_one() {
    local id=$1
    local prefix="$TREES_DIR/$id"
    [[ -s "$prefix.iqtree" && -s "$prefix.treefile" ]] && return 0
    if ! iqtree3 -s "$TRIM_DIR/$id.faa" -m MFP --msub nuclear -B 1000 -bnni -T 1 \
        --prefix "$prefix" --quiet > /dev/null 2>&1; then
        echo "  FAILED iqtree3: $id (see $prefix.log)" >&2
        return 1
    fi
    rm -f "$prefix".{ckp.gz,mldist,bionj,model.gz,uniqueseq.phy,splits.nex,contree,ufboot}
}
export -f tree_one
export TRIM_DIR TREES_DIR

grep . "$MARKERS_LIST" | xargs -P "$SLURM_CPUS_PER_TASK" -I {} bash -c 'tree_one "$1"' _ {}

N_MARKERS=$(grep -c . "$MARKERS_LIST")
N_DONE=$(grep . "$MARKERS_LIST" | while read -r id; do [[ -s "$TREES_DIR/$id.treefile" ]] && echo; done | wc -l)
echo "Finished $(date): $N_DONE of $N_MARKERS marker trees"
if (( N_DONE < N_MARKERS )); then
    echo "$(( N_MARKERS - N_DONE )) markers have no tree yet, resubmit to continue them" >&2
    exit 1
fi
