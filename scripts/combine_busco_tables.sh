#!/bin/bash
# combine_busco_tables.sh
# Description: Concatenates the BUSCO full_table.tsv of every proteome into a single table,
#              with 'portal' and 'batch' columns in front and BUSCO's comment lines removed, so
#              all results can be filtered and compared at once. Choosing the batch runs keeps
#              separate tables for different proteome versions, e.g. before and after collapsing
#              isoforms.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-14
# Version: 1.0
#
# Usage (from the project folder):
#   scripts/combine_busco_tables.sh [batches] [portal_list] [results_dir] [output_file]
#   batches      Comma-separated batch runs (folders in results_dir) to combine, or 'all' for every
#                batch run (default: all). For example:
#                  scripts/combine_busco_tables.sh proteomes,rerun        (before collapsing isoforms)
#                  scripts/combine_busco_tables.sh collapsed_proteomes    (after collapsing isoforms)
#   portal_list  TSV with a header row and portal IDs in column 1
#                (default: 00_inputs/mycocosm_datasets.tsv).
#   results_dir  BUSCO results folder (default: 04_BUSCO/results).
#   output_file  Where the combined table is saved
#                (default: 04_BUSCO/busco_full_table_<batches>.tsv, e.g. busco_full_table_proteomes_rerun.tsv).
#
# Expects BUSCO batch-mode output: <results_dir>/<batch>/<portal>_proteins.faa/run_<lineage>/full_table.tsv
# Each portal must be found in exactly one of the chosen batches.

BATCHES=${1:-all}
PORTAL_LIST=${2:-00_inputs/mycocosm_datasets.tsv}
RESULTS_DIR=${3:-04_BUSCO/results}
OUT_FILE=${4:-04_BUSCO/busco_full_table_${BATCHES//,/_}.tsv}

if [[ ! -f "$PORTAL_LIST" ]]; then
    echo "Portal list file '$PORTAL_LIST' not found" >&2
    echo "Usage: $0 [batches] [portal_list] [results_dir] [output_file]" >&2
    exit 1
fi

if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "BUSCO results folder '$RESULTS_DIR' not found" >&2
    echo "Usage: $0 [batches] [portal_list] [results_dir] [output_file]" >&2
    exit 1
fi

# Batch runs to search: every folder in results_dir, or the comma-separated list
if [[ "$BATCHES" == all ]]; then
    BATCH_DIRS=("$RESULTS_DIR"/*/)
else
    IFS=, read -ra BATCH_NAMES <<< "$BATCHES"
    BATCH_DIRS=()
    for BATCH in "${BATCH_NAMES[@]}"; do
        if [[ ! -d "$RESULTS_DIR/$BATCH" ]]; then
            echo "Batch run '$RESULTS_DIR/$BATCH' not found; batch runs there:" >&2
            find "$RESULTS_DIR" -mindepth 1 -maxdepth 1 -type d -printf '  %f\n' >&2
            exit 1
        fi
        BATCH_DIRS+=("$RESULTS_DIR/$BATCH/")
    done
fi

# Convert Windows line endings in place, otherwise every portal ID ends in a hidden \r
sed -i 's/\r$//' "$PORTAL_LIST"

# Portal IDs: column 1, skipping the header row
PORTALS=$(awk -F '\t' 'NR > 1 && $1 != "" { print $1 }' "$PORTAL_LIST")

TABLES=()
MISSING=()
DUPLICATED=()
for PORTAL in $PORTALS; do
    FOUND=()
    for BATCH_DIR in "${BATCH_DIRS[@]}"; do
        for TABLE in "$BATCH_DIR${PORTAL}_proteins.faa"/run_*/full_table.tsv; do
            [[ -f "$TABLE" ]] && FOUND+=("$TABLE")
        done
    done
    (( ${#FOUND[@]} )) || { MISSING+=("$PORTAL"); continue; }
    (( ${#FOUND[@]} == 1 )) || DUPLICATED+=("$PORTAL: ${FOUND[*]}")
    TABLES+=("${FOUND[0]}")
done

# A portal with results in two of the chosen batch runs is ambiguous, so the table is not written
if (( ${#DUPLICATED[@]} )); then
    echo "${#DUPLICATED[@]} portals have a full_table.tsv in more than one batch run;" \
        "choose the batch runs to combine (see Usage):" >&2
    printf '  %s\n' "${DUPLICATED[@]}" >&2
    exit 1
fi

if (( ${#TABLES[@]} == 0 )); then
    echo "No full_table.tsv found in batch runs '$BATCHES' for the portals in '$PORTAL_LIST'" >&2
    exit 1
fi

# BUSCO's own column header is the comment line starting with '# Busco id'; it is written once,
# with 'portal' and 'batch' in front. Other comment lines (BUSCO version, lineage) are dropped.
# Writes to a .part file first, so a failed run never overwrites an existing table.
mkdir -p "$(dirname "$OUT_FILE")"
N_DONE=0
for TABLE in "${TABLES[@]}"; do
    PORTAL=$(basename "$(dirname "$(dirname "$TABLE")")" _proteins.faa)
    BATCH=$(basename "$(dirname "$(dirname "$(dirname "$TABLE")")")")
    N_DONE=$(( N_DONE + 1 ))
    awk -F '\t' -v OFS='\t' -v portal="$PORTAL" -v batch="$BATCH" -v first=$(( N_DONE == 1 )) '
        /^# Busco id/ { if (first) { sub(/^# /, ""); print "portal", "batch", $0 }; next }
        /^#/ || NF == 0 { next }
        { print portal, batch, $0 }' "$TABLE"
done > "$OUT_FILE.part"
mv "$OUT_FILE.part" "$OUT_FILE"

echo "Combined ${#TABLES[@]} BUSCO tables into $OUT_FILE ($(( $(wc -l < "$OUT_FILE") - 1 )) rows)"
if (( ${#MISSING[@]} )); then
    echo "${#MISSING[@]} portals on the list have no BUSCO full_table.tsv:" >&2
    printf '  %s\n' "${MISSING[@]}" >&2
fi
