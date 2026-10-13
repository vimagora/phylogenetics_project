#!/bin/bash
# list_mycocosm_genomes.sh
# Description: Uses the mycocosm catalog API to retrieve the list of all published genomes
#              (name, portal ID, NCBI taxon, assembly stats and publication) as a TSV.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-13
# Version: 1.0
#
# Usage (from the project folder): scripts/list_mycocosm_genomes.sh [output_file] [group]
#   output_file  Where the list is saved (default: 00_inputs/mycocosm_published.tsv).
#   group        mycocosm group to list (default: fungi).
#
# The catalog does not need an API key (downloading the genomes does).

OUT_FILE=${1:-00_inputs/mycocosm_published.tsv}
GROUP=${2:-fungi}
API='https://mycocosm.jgi.doe.gov/ext-api/mycocosm/catalog/download-group-tsv'

# Downloads to a .part file that is only renamed once it looks like the catalog,
# so a failed request never overwrites an existing list
PART="$OUT_FILE.part"
mkdir -p "$(dirname "$OUT_FILE")"
echo "Listing published genomes for group '$GROUP' ..."
code=$(curl --silent --get "$API?flt=&seq=all&pub=published&grp=$GROUP&srt=name" --output "$PART" --write-out '%{http_code}')
if [[ "$code" != 200 ]]; then
    echo "FAILED (HTTP $code)" >&2
    [[ -f "$PART" ]] && echo "server said: $(head -c 300 "$PART")" >&2
    rm -f "$PART"
    exit 1
fi
if ! head -n 1 "$PART" | grep -q $'\tportal\t'; then
    echo "FAILED (response is not the catalog TSV)" >&2
    echo "server said: $(head -c 300 "$PART")" >&2
    rm -f "$PART"
    exit 1
fi
mv "$PART" "$OUT_FILE"

echo "Saved $(( $(wc -l < "$OUT_FILE") - 1 )) genomes to $OUT_FILE"
