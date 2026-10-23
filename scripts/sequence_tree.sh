#!/bin/bash
# sequence_tree.sh
# Family-specific: sugar transporters (structure concept test; see README, "Sugar transporter structures").
# Description: Sequence-only baseline for the structure trees: a maximum-likelihood tree of the same
#              sugar transporters (those with a trimmed model), from their full sequences in
#              08_structures/structure_ids.tsv. MAFFT L-INS-i alignment, ClipKIT gappy trimming
#              (columns with > 50% gaps removed, as for the species tree) and IQ-TREE 3 with
#              ModelFinder and 1000 ultrafast bootstraps.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-22
# Version: 1.0
#
# Usage (from the project folder; needs mafft, clipkit and iqtree3, e.g. the structure_test env):
#   scripts/sequence_tree.sh [models_dir] [out_dir] [threads]
#   models_dir  Trimmed models; their names select the proteins (default: 08_structures/models_trimmed).
#   out_dir     Where results are saved (default: 08_structures/trees).
#   threads     Threads for MAFFT and IQ-TREE (default: all).
#
# Output, inside out_dir: sequences.faa, sequences_aligned.faa, sequences_trimmed.faa and
# sequence_ml.* (IQ-TREE; sequence_ml.treefile is the tree with UFBoot support).

MODELS_DIR=${1:-08_structures/models_trimmed}
OUT_DIR=${2:-08_structures/trees}
THREADS=${3:-$(nproc)}
IDS_TABLE=08_structures/structure_ids.tsv

for TOOL in mafft clipkit iqtree3; do
    if ! command -v "$TOOL" > /dev/null; then
        echo "$TOOL not found; activate the structure_test environment" >&2
        exit 1
    fi
done
if ! compgen -G "$MODELS_DIR/*.pdb" > /dev/null; then
    echo "No models in '$MODELS_DIR', run scripts/fetch_afdb_models.py first" >&2
    exit 1
fi

set -o pipefail
mkdir -p "$OUT_DIR"

# Sequences of the proteins that have a model
find "$MODELS_DIR" -maxdepth 1 -name '*.pdb' -printf '%f\n' | sed 's/\.pdb$//' | sort > "$OUT_DIR/ids.txt"
awk -F'\t' 'NR == FNR { keep[$1] = 1; next }
    FNR > 1 && ($1 in keep) { print ">" $1; print $9 }' "$OUT_DIR/ids.txt" "$IDS_TABLE" > "$OUT_DIR/sequences.faa"
N_SEQ=$(grep -c '^>' "$OUT_DIR/sequences.faa")
if (( N_SEQ != $(wc -l < "$OUT_DIR/ids.txt") )); then
    echo "Some models have no sequence in $IDS_TABLE" >&2
    exit 1
fi
echo "Aligning $N_SEQ sequences ..."

mafft --localpair --maxiterate 1000 --thread "$THREADS" --quiet "$OUT_DIR/sequences.faa" \
    > "$OUT_DIR/sequences_aligned.faa" || exit 1
clipkit "$OUT_DIR/sequences_aligned.faa" -m gappy -g 0.5 -of fasta -o "$OUT_DIR/sequences_trimmed.faa" \
    > /dev/null || exit 1
iqtree3 -s "$OUT_DIR/sequences_trimmed.faa" -m MFP -B 1000 -bnni -T "$THREADS" \
    --prefix "$OUT_DIR/sequence_ml" --redo --quiet || exit 1

echo "Sequence tree: $OUT_DIR/sequence_ml.treefile" \
    "($(grep -m1 'Best-fit model' "$OUT_DIR/sequence_ml.iqtree" | sed 's/.*: //'))"
