#!/bin/bash
#SBATCH --job-name=family_tree
#SBATCH --output=09_families/logs/tree_%j.out
#SBATCH --error=09_families/logs/tree_%j.err
#SBATCH --open-mode=append
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
# family_tree.sbatch
# Description: Slurm job that infers a gene family's tree from the trimmed alignment
#              (align_family.sbatch), in one of three modes, each reviewed before the next:
#              fasttree     FastTree (LG + Gamma), a quick tree to check the outgroup and the
#                           expectations (families/<family>/expectations.md) before the long runs;
#              modelfinder  IQ-TREE 3 ModelFinder (-m MF) on the full alignment;
#              tree         IQ-TREE 3 with the ModelFinder model, 1000 ultrafast bootstraps with NNI
#                           refinement and 1000 SH-aLRT replicates, drawn rooted on the outgroups
#                           marked root = TRUE in outgroups.tsv.
#              IQ-TREE picks its number of threads (-T AUTO, at most the CPUs of the job): family
#              alignments are short after trimming, and more threads than columns can use only add cost.
#              The IQ-TREE modes resume from their checkpoint when resubmitted (e.g. after the time
#              limit); a finished run is not repeated.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-23
# Version: 1.0
#
# Usage (from the project folder): sbatch [--time=...] scripts/family_tree.sbatch <family> <mode>
#   e.g. sbatch --time=04:00:00 scripts/family_tree.sbatch sugar_transporters fasttree
#
# Output, inside <OUT_DIR>/tree (OUT_DIR from families/<family>/family.conf):
#   fasttree.nwk, fasttree.log     FastTree tree (local support values) and log
#   modelfinder.*                  ModelFinder results; the best model is in modelfinder.iqtree
#   ml.*                           IQ-TREE results; ml.treefile is the tree with SH-aLRT/UFBoot support
#
# The Slurm account (and ORTHOFINDER_DIR, a fallback source of FastTree) come from config.sh:
# run 'source config.sh' before sbatch.

# Slurm starts the job in the folder it was submitted from, which must be the project folder
if [[ ! -d scripts || ! -d 00_inputs ]]; then
    echo "Submit this job from the project folder" >&2
    exit 1
fi

FAMILY=$1
MODE=$2
CONF="families/$FAMILY/family.conf"
if [[ -z "$FAMILY" || ! -s "$CONF" || ! "$MODE" =~ ^(fasttree|modelfinder|tree)$ ]]; then
    echo "Usage: sbatch scripts/family_tree.sbatch <family> fasttree|modelfinder|tree" >&2
    exit 1
fi
# shellcheck source=/dev/null
source "$CONF"
DIR="$OUT_DIR/tree"
ALIGNMENT="$DIR/trimmed.faa"
THREADS=${SLURM_CPUS_PER_TASK:-4}
if [[ ! -s "$ALIGNMENT" ]]; then
    echo "'$ALIGNMENT' not found, run scripts/align_family.sbatch $FAMILY first" >&2
    exit 1
fi

echo "Job $SLURM_JOB_ID on $(hostname), started $(date): $MODE for $NAME"
module load bio-apps/v202603

if [[ $MODE == fasttree ]]; then
    if [[ -s "$DIR/fasttree.nwk" ]]; then
        echo "FastTree tree already exists: $DIR/fasttree.nwk"
        exit 0
    fi
    # FastTree from a module, else from the OrthoFinder container (which ships it)
    module load fasttree 2> /dev/null
    if ! command -v FastTreeMP > /dev/null && ! command -v FastTree > /dev/null && ! command -v fasttree > /dev/null; then
        [[ -f config.sh ]] && source config.sh
        [[ -n "$ORTHOFINDER_DIR" ]] && export PATH="$ORTHOFINDER_DIR/bin:$PATH"
    fi
    FASTTREE=$(command -v FastTreeMP || command -v FastTree || command -v fasttree)
    if [[ -z "$FASTTREE" ]]; then
        echo "FastTree not found (no module, and not in \$ORTHOFINDER_DIR/bin)" >&2
        exit 1
    fi
    export OMP_NUM_THREADS=$THREADS
    if ! "$FASTTREE" -lg -gamma -log "$DIR/fasttree.log" "$ALIGNMENT" > "$DIR/fasttree.nwk.part"; then
        echo "FastTree failed, see $DIR/fasttree.log" >&2
        rm -f "$DIR/fasttree.nwk.part"
        exit 1
    fi
    mv "$DIR/fasttree.nwk.part" "$DIR/fasttree.nwk"
    echo "Done $(date): $DIR/fasttree.nwk ($FASTTREE). Next: python3 scripts/make_family_itol.py $FAMILY, then review"
    exit 0
fi

module load iqtree3/3.1.4
if ! command -v iqtree3 > /dev/null; then
    echo "iqtree3 not found, check that the modules above loaded" >&2
    exit 1
fi

if [[ $MODE == modelfinder ]]; then
    PREFIX="$DIR/modelfinder"
    if [[ -s "$PREFIX.iqtree" ]]; then
        echo "ModelFinder already finished: $(grep -m1 'Best-fit model' "$PREFIX.iqtree")"
        exit 0
    fi
    iqtree3 -s "$ALIGNMENT" -m MF -T AUTO --threads-max "$THREADS" --prefix "$PREFIX" --quiet
    STATUS=$?
    (( STATUS == 0 )) && echo "Done $(date): $(grep -m1 'Best-fit model' "$PREFIX.iqtree")"
    exit $STATUS
fi

# tree: the ModelFinder model, rooted on the outgroups marked for rooting
MODEL=$(grep -m1 'Best-fit model according to BIC' "$DIR/modelfinder.iqtree" 2> /dev/null | sed 's/.*: //')
if [[ -z "$MODEL" ]]; then
    echo "No ModelFinder result in $DIR/modelfinder.iqtree; run the modelfinder mode first" >&2
    exit 1
fi
ROOT=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
    $col["root"] == "TRUE" { printf "%s%s", sep, $col["outgroup_id"]; sep = "," }' "families/$FAMILY/outgroups.tsv")
PREFIX="$DIR/ml"
if [[ -s "$PREFIX.iqtree" && -s "$PREFIX.treefile" ]]; then
    echo "Tree already finished: $PREFIX.treefile"
    exit 0
fi
echo "Model $MODEL; rooted on $ROOT"
iqtree3 -s "$ALIGNMENT" -m "$MODEL" -B 1000 -bnni -alrt 1000 ${ROOT:+-o "$ROOT"} -T AUTO --threads-max "$THREADS" \
    --prefix "$PREFIX" --quiet
STATUS=$?
echo "Finished $(date) with exit code $STATUS"
(( STATUS == 0 )) && echo "Tree: $PREFIX.treefile. Next: python3 scripts/make_family_itol.py $FAMILY, then review"
exit $STATUS
