#!/bin/bash
#SBATCH --job-name=species_trees
#SBATCH --output=06_species_tree/logs/species_trees_%j.out
#SBATCH --error=06_species_tree/logs/species_trees_%j.err
#SBATCH --open-mode=append
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
# species_trees.sbatch
# Description: Slurm job that infers the two species trees from the kept BUSCO markers:
#              1. ASTRAL-IV (ASTER) on the per-marker trees: the coalescent species tree, with
#                 local posterior probabilities; takes minutes;
#              2. IQ-TREE 3 on the supermatrix, one partition per marker with the marker's best-fit
#                 model (edge-linked proportional branch lengths), 1000 ultrafast bootstraps with
#                 NNI refinement and 1000 SH-aLRT replicates; takes hours to days.
#              Both trees are rooted on the outgroup; ASTER recommends correct rooting for accurate
#              ASTRAL branch lengths. The default, Catan2 (Catenaria anguillulae), puts
#              Blastocladiomycota first, as in the concatenated tree of Li et al. 2021; the cohort
#              has no non-fungal genome. (The first runs used Spipu1, which falls inside
#              Chytridiomyceta and makes it look paraphyletic.)
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-20
# Version: 1.1
#
# Usage (from the project folder): sbatch scripts/species_trees.sbatch [--outgroup PORTAL]
#   --outgroup  One cohort portal to root both trees on (default: Catan2).
# To rerun only ASTRAL, e.g. with another outgroup, move the old astral/ folder away and resubmit
# with a short time limit: the finished IQ-TREE run is skipped.
#   mv 06_species_tree/astral 06_species_tree/astral_old
#   sbatch --time=01:00:00 scripts/species_trees.sbatch --outgroup Catan2
# Inputs: 06_species_tree/marker_trees.nwk and concatenated/partitions_models.nex from
# collect_marker_trees.py, and concatenated/supermatrix.faa from concatenate_markers.py.
#
# Output, inside 06_species_tree:
#   astral/species_tree_astral.nwk   ASTRAL-IV tree (internal branch lengths in coalescent units)
#   concatenated/iqtree/supermatrix.*  IQ-TREE results; .treefile is the tree with SH-aLRT/UFBoot support
# A finished ASTRAL tree is not recomputed. IQ-TREE resumes from its checkpoint if the job is
# resubmitted (e.g. after the time limit), so resubmitting continues the run.
#
# The Slurm account comes from config.sh: run 'source config.sh' before sbatch.

# Slurm starts the job in the folder it was submitted from, which must be the project folder
if [[ ! -d scripts || ! -d 00_inputs ]]; then
    echo "Submit this job from the project folder" >&2
    exit 1
fi

DIR=06_species_tree
GENE_TREES="$DIR/marker_trees.nwk"
SUPERMATRIX="$DIR/concatenated/supermatrix.faa"
PARTITIONS="$DIR/concatenated/partitions_models.nex"
ASTRAL_TREE="$DIR/astral/species_tree_astral.nwk"
IQTREE_PREFIX="$DIR/concatenated/iqtree/supermatrix"
OUTGROUP=Catan2

while (( $# > 0 )); do
    case "$1" in
        --outgroup)
            if [[ -z "${2:-}" ]]; then echo "--outgroup needs a portal" >&2; exit 1; fi
            OUTGROUP=$2; shift 2 ;;
        *) echo "Unknown argument '$1'; usage: sbatch scripts/species_trees.sbatch [--outgroup PORTAL]" >&2; exit 1 ;;
    esac
done
if ! awk -F'\t' -v p="$OUTGROUP" 'NR > 1 { sub(/\r$/, "", $1); if ($1 == p) found = 1 } END { exit !found }' \
    00_inputs/cohort_datasets.tsv; then
    echo "Outgroup '$OUTGROUP' is not a cohort portal (00_inputs/cohort_datasets.tsv)" >&2
    exit 1
fi

for FILE in "$GENE_TREES" "$SUPERMATRIX" "$PARTITIONS"; do
    if [[ ! -s "$FILE" ]]; then
        echo "'$FILE' not found, run scripts/concatenate_markers.py and scripts/collect_marker_trees.py first" >&2
        exit 1
    fi
done

module load bio-apps/v202603
module load aster/1.25
module load iqtree3/3.1.4
# ASTRAL-IV is 'astral4' in recent ASTER releases, 'astral' in older ones
ASTRAL=$(command -v astral4 || command -v astral)
if [[ -z "$ASTRAL" ]] || ! command -v iqtree3 > /dev/null; then
    echo "astral4/astral or iqtree3 not found, check that the modules above loaded" >&2
    exit 1
fi

echo "Job $SLURM_JOB_ID on $(hostname), started $(date); outgroup $OUTGROUP"

# 1. ASTRAL-IV on the per-marker trees
if [[ -s "$ASTRAL_TREE" ]]; then
    echo "ASTRAL tree already exists: $ASTRAL_TREE"
else
    echo "Running $(basename "$ASTRAL") on $(grep -c . "$GENE_TREES") marker trees ..."
    mkdir -p "$(dirname "$ASTRAL_TREE")"
    if ! "$ASTRAL" -i "$GENE_TREES" -o "$ASTRAL_TREE.part" --root "$OUTGROUP" -t "$SLURM_CPUS_PER_TASK" \
        > "$DIR/astral/astral.out" 2>&1; then
        echo "ASTRAL failed, see $DIR/astral/astral.out" >&2
        rm -f "$ASTRAL_TREE.part"
        exit 1
    fi
    mv "$ASTRAL_TREE.part" "$ASTRAL_TREE"
    echo "ASTRAL tree: $ASTRAL_TREE"
fi

# 2. IQ-TREE 3 on the supermatrix, resuming from the checkpoint if there is one. A finished run
#    (with its .iqtree report) is not repeated, since IQ-TREE refuses to overwrite it.
if [[ -s "$IQTREE_PREFIX.iqtree" && -s "$IQTREE_PREFIX.treefile" ]]; then
    echo "Concatenated tree already exists: $IQTREE_PREFIX.treefile (written with the outgroup of its own run;"
    echo "the root of an IQ-TREE tree only sets how it is drawn, so reroot it when plotting)"
    exit 0
fi
mkdir -p "$(dirname "$IQTREE_PREFIX")"
echo "Running IQ-TREE 3 on $SUPERMATRIX ($(date)) ..."
iqtree3 -s "$SUPERMATRIX" -p "$PARTITIONS" -B 1000 -bnni -alrt 1000 -o "$OUTGROUP" \
    -T "$SLURM_CPUS_PER_TASK" --prefix "$IQTREE_PREFIX" --quiet
IQTREE_STATUS=$?

echo "Finished $(date) with exit code $IQTREE_STATUS"
(( IQTREE_STATUS == 0 )) && echo "Concatenated tree: $IQTREE_PREFIX.treefile"
exit $IQTREE_STATUS
