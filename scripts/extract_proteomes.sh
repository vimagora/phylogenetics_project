#!/bin/bash
# extract_proteomes.sh
# Description: Decompresses all downloaded proteome files (<portal>_proteins.faa.gz)
#              into a single folder, as input for BUSCO.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-13
# Version: 1.0
#
# Usage (from the project folder): scripts/extract_proteomes.sh [portal_list] [genomes_dir] [output_dir]
#   portal_list  TSV with a header row and portal IDs in column 1
#                (default: 00_inputs/mycocosm_datasets.tsv).
#   genomes_dir  Folder with one subdirectory per portal (default: 01_genomes).
#   output_dir   Where the decompressed proteomes are saved (default: 02_BUSCO/proteomes).
#
# Proteomes that were already decompressed are skipped, so it is safe to re-run
# after downloading more genomes.

PORTAL_LIST=${1:-00_inputs/mycocosm_datasets.tsv}
GENOMES_DIR=${2:-01_genomes}
OUT_DIR=${3:-02_BUSCO/proteomes}

if [[ ! -f "$PORTAL_LIST" ]]; then
    echo "Portal list file '$PORTAL_LIST' not found" >&2
    echo "Usage: $0 [portal_list] [genomes_dir] [output_dir]" >&2
    exit 1
fi

if [[ ! -d "$GENOMES_DIR" ]]; then
    echo "Genomes folder '$GENOMES_DIR' not found" >&2
    echo "Usage: $0 [portal_list] [genomes_dir] [output_dir]" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

# Convert Windows line endings in place, otherwise every portal ID ends in a hidden \r
sed -i 's/\r$//' "$PORTAL_LIST"

# Portal IDs: column 1, skipping the header row
PORTALS=$(awk -F '\t' 'NR > 1 && $1 != "" { print $1 }' "$PORTAL_LIST")

set -o pipefail

N_EXTRACTED=0
N_SKIPPED=0
FAILED=()
MISSING=()
for PORTAL in $PORTALS; do
    PROTEOME="$GENOMES_DIR"/"$PORTAL"/"$PORTAL"_proteins.faa.gz

    # Portals on the list without a proteome (e.g. a failed or unfinished download)
    [[ -f "$PROTEOME" ]] || MISSING+=("$PORTAL")
    [[ -f "$PROTEOME" ]] || continue

    OUT="$OUT_DIR/$(basename "$PROTEOME" .gz)"
    if [[ -s "$OUT" ]]; then
        N_SKIPPED=$(( N_SKIPPED + 1 ))
        continue
    fi
    # Decompress to a .part file first, so an interrupted run never leaves a truncated proteome
    if gzip -dc "$PROTEOME" | sed '/^;/d' > "$OUT.part"; then
        mv "$OUT.part" "$OUT"
        N_EXTRACTED=$(( N_EXTRACTED + 1 ))
    else
        rm -f "$OUT.part"
        FAILED+=("$PROTEOME")
    fi
done

echo "Decompressed $N_EXTRACTED proteomes into $OUT_DIR ($N_SKIPPED already there)"
if (( ${#MISSING[@]} )); then
    echo "${#MISSING[@]} portals on the list have no proteome:" >&2
    printf '  %s\n' "${MISSING[@]}" >&2
fi
if (( ${#FAILED[@]} )); then
    echo "${#FAILED[@]} proteomes could not be decompressed:" >&2
    printf '  %s\n' "${FAILED[@]}" >&2
    exit 1
fi
