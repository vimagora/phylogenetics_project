#!/bin/bash
# download_mycocosm_genomes.sh
# Description: Uses the mycocosm API to download the gene models (GFF3), repeat-masked
#              assembly and default protein and CDS sequences (FASTA) for a set of portals.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-13
# Version: 1.0
#
# Usage (from the project folder): scripts/download_mycocosm_genomes.sh [portal_list] [output_dir]
#   portal_list  TSV with a header row and portal IDs in column 1
#                (default: 00_inputs/mycocosm_datasets.tsv).
#   output_dir   Where files are saved (default: 01_genomes).
#                Each portal gets its own subdirectory.
#
# Requires a JGI API key, saved alone in ~/.jgi_api_key (or set MYCOCOSM_API_KEY_FILE).
# It is sent as the 'Authorization' header; without it every download returns HTTP 403.

PORTAL_LIST=${1:-00_inputs/mycocosm_datasets.tsv}
OUT_DIR=${2:-01_genomes}
API='https://mycocosm.jgi.doe.gov/ext-api/mycocosm/download'
API_KEY_FILE=${MYCOCOSM_API_KEY_FILE:-$HOME/.jgi_api_key}

if [[ ! -f "$PORTAL_LIST" ]]; then
    echo "Portal list '$PORTAL_LIST' not found" >&2
    echo "Usage: $0 [portal_list] [output_dir]" >&2
    exit 1
fi

if [[ ! -s "$API_KEY_FILE" ]]; then
    echo "API key file '$API_KEY_FILE' not found or empty" >&2
    exit 1
fi
API_KEY=$(<"$API_KEY_FILE")

# Downloads one file, skipping it if it already exists. Downloads go to a .part file that is
# only renamed once it is a complete gzip, so an interrupted run never leaves a truncated file.
# Stops the whole script if the key is refused, since every other download would fail too.
download() {
    local url=$1 out=$2 part="$2.part" code
    if [[ -s "$out" ]]; then
        echo "  skipping $(basename "$out") (already exists)"
        return 0
    fi
    echo "  downloading $(basename "$out") ..."
    code=$(curl --silent --get "$url" -H "Authorization: $API_KEY" --output "$part" --write-out '%{http_code}')
    if [[ "$code" != 200 ]]; then
        echo "  FAILED (HTTP $code): $url" >&2
        [[ -f "$part" ]] && echo "  server said: $(head -c 300 "$part")" >&2
        rm -f "$part"
        if [[ "$code" == 401 || "$code" == 403 ]]; then
            echo "API key refused or expired ($API_KEY_FILE), stopping" >&2
            exit 1
        fi
        return 1
    fi
    if ! gzip -t "$part" 2>/dev/null; then
        echo "  FAILED (not a complete gzip file): $url" >&2
        grep -Iq . "$part" && echo "  server said: $(head -c 300 "$part")" >&2
        rm -f "$part"
        return 1
    fi
    mv "$part" "$out"
}

# Convert Windows line endings in place, otherwise every portal ID ends in a hidden \r
dos2unix --quiet "$PORTAL_LIST"

# Portal IDs: column 1, skipping the header row
PORTALS=$(awk -F '\t' 'NR > 1 && $1 != "" { print $1 }' "$PORTAL_LIST")

FAILED=()
for PORTAL in $PORTALS; do
    echo "Portal '$PORTAL'"
    mkdir -p "$OUT_DIR/$PORTAL"
    download "$API/$PORTAL/genes-gff3" "$OUT_DIR/$PORTAL/${PORTAL}_genes.gff3.gz" || FAILED+=("$PORTAL genes")
    download "$API/$PORTAL/repeatmasked-assembly" "$OUT_DIR/$PORTAL/${PORTAL}_repeatmasked-assembly.fasta.gz" || FAILED+=("$PORTAL assembly")
    download "$API/$PORTAL/protein-default" "$OUT_DIR/$PORTAL/${PORTAL}_proteins.faa.gz" || FAILED+=("$PORTAL proteins")
    download "$API/$PORTAL/cds-default" "$OUT_DIR/$PORTAL/${PORTAL}_cds.fna.gz" || FAILED+=("$PORTAL cds")
done

if (( ${#FAILED[@]} )); then
    echo "Done with ${#FAILED[@]} failed downloads:" >&2
    printf '  %s\n' "${FAILED[@]}" >&2
    exit 1
fi
echo "All done"
