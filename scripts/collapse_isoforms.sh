#!/bin/bash
# collapse_isoforms.sh
# Description: Collapses isoforms and redundant gene models (step A of the proteome cleanup)
#              for every portal on the list, with collapse_isoforms.py. MycoCosm GFFs store
#              every isoform as a separate gene, so models that share most of their coding
#              sequence (CDS) on the same strand are grouped into one locus, and only the longest
#              protein of each locus is kept. See collapse_isoforms.py for the exact rules.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-15
# Version: 1.0
#
# Usage (from the project folder):
#   scripts/collapse_isoforms.sh [portal_list] [genomes_dir] [proteomes_dir] [output_dir]
#   portal_list    TSV with a header row and portal IDs in column 1
#                  (default: 00_inputs/mycocosm_datasets.tsv).
#   genomes_dir    Folder with one subdirectory per portal, with <portal>_genes.gff3.gz
#                  (default: 01_genomes).
#   proteomes_dir  Decompressed proteomes from extract_proteomes.sh (default: 03_BUSCO/proteomes).
#   output_dir     Where the results are saved (default: 04_isoforms).
#
# Output, inside output_dir:
#   collapsed_proteomes/<portal>_proteins.faa  one protein per locus
#   removed_models.tsv                         every removed model and the model kept in its place
#   summary.tsv                                models before and after, per portal
#   unmerged_overlaps.tsv                      overlapping pairs that were NOT merged, for review
#   work/<portal>/                             intermediate files, kept for checking
#
# Portals that already have a collapsed proteome are skipped.
# Uses the python-data module on Roihu (pandas, numpy, networkx).

PORTAL_LIST=${1:-00_inputs/mycocosm_datasets.tsv}
GENOMES_DIR=${2:-01_genomes}
PROTEOMES_DIR=${3:-03_BUSCO/proteomes}
OUT_DIR=${4:-04_isoforms}
COLLAPSE_PY="$(dirname "$0")/collapse_isoforms.py"

# Minimum shared CDS, as % of the shorter model's CDS, for two models to be merged. On Neucr2,
# every Broad isoform pair shares >= 50% (80% would split 9 real loci); on Agabi_varbisH97_2,
# end-overlapping neighbours share <= 5% and one partial overlap shares 33%.
MIN_SHARED_PCT=50

if [[ ! -f "$PORTAL_LIST" ]]; then
    echo "Portal list file '$PORTAL_LIST' not found" >&2
    echo "Usage: $0 [portal_list] [genomes_dir] [proteomes_dir] [output_dir]" >&2
    exit 1
fi

for DIR in "$GENOMES_DIR" "$PROTEOMES_DIR"; do
    if [[ ! -d "$DIR" ]]; then
        echo "Folder '$DIR' not found" >&2
        echo "Usage: $0 [portal_list] [genomes_dir] [proteomes_dir] [output_dir]" >&2
        exit 1
    fi
done

module load python-data
if ! python3 -c "import pandas, numpy, networkx" 2> /dev/null; then
    echo "python3 with pandas, numpy and networkx not found, check that python-data loaded" >&2
    exit 1
fi

# Convert Windows line endings in place, otherwise every portal ID ends in a hidden \r
sed -i 's/\r$//' "$PORTAL_LIST"

# Portal IDs: column 1, skipping the header row
PORTALS=$(awk -F '\t' 'NR > 1 && $1 != "" { print $1 }' "$PORTAL_LIST")

mkdir -p "$OUT_DIR/collapsed_proteomes" "$OUT_DIR/work"

N_DONE=0
N_SKIPPED=0
FAILED=()
MISSING=()
for PORTAL in $PORTALS; do
    GFF="$GENOMES_DIR/$PORTAL/${PORTAL}_genes.gff3.gz"
    PROTEOME="$PROTEOMES_DIR/${PORTAL}_proteins.faa"
    OUT="$OUT_DIR/collapsed_proteomes/${PORTAL}_proteins.faa"
    if [[ ! -f "$GFF" || ! -s "$PROTEOME" ]]; then
        MISSING+=("$PORTAL")
        continue
    fi
    # Skip portals already collapsed with this version of the script (summary row with 11 columns)
    if [[ -s "$OUT" ]] && awk -F '\t' 'NF == 11 { ok = 1 } END { exit !ok }' \
        "$OUT_DIR/work/$PORTAL/summary_row.tsv" 2> /dev/null; then
        N_SKIPPED=$(( N_SKIPPED + 1 ))
        continue
    fi
    if python3 "$COLLAPSE_PY" --portal "$PORTAL" --gff "$GFF" --proteome "$PROTEOME" \
        --work-dir "$OUT_DIR/work/$PORTAL" --out-fasta "$OUT" --min-shared-pct "$MIN_SHARED_PCT"; then
        N_DONE=$(( N_DONE + 1 ))
    else
        FAILED+=("$PORTAL")
    fi
done

# Combined tables for every portal on the list that has been collapsed, rebuilt on every run
{
    printf 'portal\tproteins_in\tgff_models\tgff_without_protein\tloci\tremoved\tpct_removed\tnot_in_gff\tproteins_out'
    printf '\tunmerged_pairs\tmax_unmerged_pct\n'
    for PORTAL in $PORTALS; do
        [[ -f "$OUT_DIR/work/$PORTAL/summary_row.tsv" ]] && cat "$OUT_DIR/work/$PORTAL/summary_row.tsv"
    done
} > "$OUT_DIR/summary.tsv"
{
    printf 'portal\tremoved_proteinId\tkept_proteinId\tremoved_length\tkept_length\n'
    for PORTAL in $PORTALS; do
        [[ -f "$OUT_DIR/work/$PORTAL/summary_row.tsv" ]] && cat "$OUT_DIR/work/$PORTAL/removed.tsv"
    done
} > "$OUT_DIR/removed_models.tsv"
{
    printf 'portal\tproteinId_a\tproteinId_b\tshared_cds_bp\tpct_of_shorter_cds\n'
    for PORTAL in $PORTALS; do
        [[ -f "$OUT_DIR/work/$PORTAL/summary_row.tsv" ]] && cat "$OUT_DIR/work/$PORTAL/unmerged.tsv"
    done
} > "$OUT_DIR/unmerged_overlaps.tsv"

echo "Collapsed $N_DONE proteomes into $OUT_DIR/collapsed_proteomes ($N_SKIPPED already there)"
echo "Summary: $OUT_DIR/summary.tsv, removed models: $OUT_DIR/removed_models.tsv," \
    "overlaps left unmerged: $OUT_DIR/unmerged_overlaps.tsv"
if (( ${#MISSING[@]} )); then
    echo "${#MISSING[@]} portals on the list have no GFF or no extracted proteome:" >&2
    printf '  %s\n' "${MISSING[@]}" >&2
fi
if (( ${#FAILED[@]} )); then
    echo "${#FAILED[@]} portals could not be collapsed:" >&2
    printf '  %s\n' "${FAILED[@]}" >&2
    exit 1
fi
