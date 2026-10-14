#!/bin/bash
# extract_proteomes.sh
# Description: Decompresses all downloaded proteome files (<portal>_proteins.faa.gz)
#              into a single folder, as input for BUSCO.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-13
# Version: 1.0
#
# Usage (from the project folder): scripts/extract_proteomes.sh [genomes_dir] [output_dir]
#   genomes_dir  Folder with one subdirectory per portal (default: 01_genomes).
#   output_dir   Where the decompressed proteomes are saved (default: 02_BUSCO/proteomes).
#
# Proteomes that were already decompressed are skipped, so it is safe to re-run
# after downloading more genomes.

GENOMES_DIR=${1:-01_genomes}
OUT_DIR=${2:-02_BUSCO/proteomes}

if [[ ! -d "$GENOMES_DIR" ]]; then
    echo "Genomes folder '$GENOMES_DIR' not found" >&2
    echo "Usage: $0 [genomes_dir] [output_dir]" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

N_EXTRACTED=0
N_SKIPPED=0
FAILED=()
for PROTEOME in "$GENOMES_DIR"/*/*_proteins.faa.gz; do
    [[ -f "$PROTEOME" ]] || continue
    OUT="$OUT_DIR/$(basename "$PROTEOME" .gz)"
    if [[ -s "$OUT" ]]; then
        N_SKIPPED=$(( N_SKIPPED + 1 ))
        continue
    fi
    # Decompress to a .part file first, so an interrupted run never leaves a truncated proteome
    if gzip -dc "$PROTEOME" > "$OUT.part"; then
        mv "$OUT.part" "$OUT"
        N_EXTRACTED=$(( N_EXTRACTED + 1 ))
    else
        rm -f "$OUT.part"
        FAILED+=("$PROTEOME")
    fi
done

# Portal folders without a proteome (e.g. a failed or unfinished download)
MISSING=()
for PORTAL_DIR in "$GENOMES_DIR"/*/; do
    PORTAL=$(basename "$PORTAL_DIR")
    [[ -f "$PORTAL_DIR/${PORTAL}_proteins.faa.gz" ]] || MISSING+=("$PORTAL")
done

echo "Decompressed $N_EXTRACTED proteomes into $OUT_DIR ($N_SKIPPED already there)"
if (( ${#MISSING[@]} )); then
    echo "${#MISSING[@]} portal folders have no proteome:" >&2
    printf '  %s\n' "${MISSING[@]}" >&2
fi
if (( ${#FAILED[@]} )); then
    echo "${#FAILED[@]} proteomes could not be decompressed:" >&2
    printf '  %s\n' "${FAILED[@]}" >&2
    exit 1
fi
