#!/bin/bash
# combine_busco_tables.sh
# Description: Concatenates the BUSCO full_table.tsv of every proteome into a single table,
#              with a 'portal' column in front and BUSCO's comment lines removed, so all
#              results can be filtered and compared at once.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-14
# Version: 1.0
#
# Usage (from the project folder): scripts/combine_busco_tables.sh [portal_list] [results_dir] [output_file]
#   portal_list  TSV with a header row and portal IDs in column 1
#                (default: 00_inputs/mycocosm_datasets.tsv).
#   results_dir  BUSCO results folder; every batch run inside it is included, e.g. both
#                proteomes/ and rerun/ (default: 02_BUSCO/results).
#   output_file  Where the combined table is saved (default: 02_BUSCO/busco_full_table_all.tsv).
#
# Expects BUSCO batch-mode output: <results_dir>/<batch>/<portal>_proteins.faa/run_<lineage>/full_table.tsv

PORTAL_LIST=${1:-00_inputs/mycocosm_datasets.tsv}
RESULTS_DIR=${2:-02_BUSCO/results}
OUT_FILE=${3:-02_BUSCO/busco_full_table_all.tsv}

if [[ ! -f "$PORTAL_LIST" ]]; then
    echo "Portal list file '$PORTAL_LIST' not found" >&2
    echo "Usage: $0 [portal_list] [results_dir] [output_file]" >&2
    exit 1
fi

if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "BUSCO results folder '$RESULTS_DIR' not found" >&2
    echo "Usage: $0 [portal_list] [results_dir] [output_file]" >&2
    exit 1
fi

# Convert Windows line endings in place, otherwise every portal ID ends in a hidden \r
sed -i 's/\r$//' "$PORTAL_LIST"

# Portal IDs: column 1, skipping the header row
PORTALS=$(awk -F '\t' 'NR > 1 && $1 != "" { print $1 }' "$PORTAL_LIST")

TABLES=()
MISSING=()
DUPLICATED=()
for PORTAL in $PORTALS; do
    FOUND=("$RESULTS_DIR"/*/"${PORTAL}_proteins.faa"/run_*/full_table.tsv)
    [[ -f "${FOUND[0]}" ]] || { MISSING+=("$PORTAL"); continue; }
    (( ${#FOUND[@]} == 1 )) || DUPLICATED+=("$PORTAL: ${FOUND[*]}")
    TABLES+=("${FOUND[0]}")
done

# A portal with results in two batch runs is ambiguous, so the table is not written
if (( ${#DUPLICATED[@]} )); then
    echo "${#DUPLICATED[@]} portals have more than one full_table.tsv, remove the outdated ones:" >&2
    printf '  %s\n' "${DUPLICATED[@]}" >&2
    exit 1
fi

if (( ${#TABLES[@]} == 0 )); then
    echo "No full_table.tsv found in '$RESULTS_DIR' for the portals in '$PORTAL_LIST'" >&2
    exit 1
fi

# BUSCO's own column header is the comment line starting with '# Busco id'; it is written once,
# with 'portal' in front. Other comment lines (BUSCO version, lineage) are dropped.
# Writes to a .part file first, so a failed run never overwrites an existing table.
mkdir -p "$(dirname "$OUT_FILE")"
N_DONE=0
for TABLE in "${TABLES[@]}"; do
    PORTAL=$(basename "$(dirname "$(dirname "$TABLE")")" _proteins.faa)
    N_DONE=$(( N_DONE + 1 ))
    awk -F '\t' -v OFS='\t' -v portal="$PORTAL" -v first=$(( N_DONE == 1 )) '
        /^# Busco id/ { if (first) { sub(/^# /, ""); print "portal", $0 }; next }
        /^#/ || NF == 0 { next }
        { print portal, $0 }' "$TABLE"
done > "$OUT_FILE.part"
mv "$OUT_FILE.part" "$OUT_FILE"

echo "Combined ${#TABLES[@]} BUSCO tables into $OUT_FILE ($(( $(wc -l < "$OUT_FILE") - 1 )) rows)"
if (( ${#MISSING[@]} )); then
    echo "${#MISSING[@]} portals on the list have no BUSCO full_table.tsv:" >&2
    printf '  %s\n' "${MISSING[@]}" >&2
fi
