#!/bin/bash
# assembly_stats.sh
# Description: Uses seqkit to make a report of the assembly statistics (number of scaffolds,
#              total length, N50, L50, GC%, number of Ns, ...) of the downloaded repeat-masked
#              assemblies. The assemblies are soft-masked, so masking does not affect the stats.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-14
# Version: 1.0
#
# Usage (from the project folder): scripts/assembly_stats.sh [portal_list] [genomes_dir] [output_file]
#   portal_list  TSV with a header row and portal IDs in column 1
#                (default: 00_inputs/mycocosm_datasets.tsv).
#   genomes_dir  Folder with one subdirectory per portal (default: 01_genomes).
#   output_file  Where the report is saved (default: 03_assembly_stats/assembly_stats.tsv).
#
# Uses SeqKit from the Roihu bio-apps modules. Uses $SLURM_CPUS_PER_TASK threads inside a
# Slurm job, 4 otherwise.

PORTAL_LIST=${1:-00_inputs/mycocosm_datasets.tsv}
GENOMES_DIR=${2:-01_genomes}
OUT_FILE=${3:-03_assembly_stats/assembly_stats.tsv}
THREADS=${SLURM_CPUS_PER_TASK:-4}

if [[ ! -f "$PORTAL_LIST" ]]; then
    echo "Portal list file '$PORTAL_LIST' not found" >&2
    echo "Usage: $0 [portal_list] [genomes_dir] [output_file]" >&2
    exit 1
fi

if [[ ! -d "$GENOMES_DIR" ]]; then
    echo "Genomes folder '$GENOMES_DIR' not found" >&2
    echo "Usage: $0 [portal_list] [genomes_dir] [output_file]" >&2
    exit 1
fi

module load bio-apps/v202603
module load seqkit/2.13.0
if ! command -v seqkit > /dev/null; then
    echo "seqkit not found, check that the modules above loaded" >&2
    exit 1
fi

# Convert Windows line endings in place, otherwise every portal ID ends in a hidden \r
sed -i 's/\r$//' "$PORTAL_LIST"

# Portal IDs: column 1, skipping the header row
PORTALS=$(awk -F '\t' 'NR > 1 && $1 != "" { print $1 }' "$PORTAL_LIST")

ASSEMBLIES=()
MISSING=()
for PORTAL in $PORTALS; do
    ASSEMBLY="$GENOMES_DIR/$PORTAL/${PORTAL}_repeatmasked-assembly.fasta.gz"
    if [[ -f "$ASSEMBLY" ]]; then
        ASSEMBLIES+=("$ASSEMBLY")
    else
        MISSING+=("$PORTAL")
    fi
done

if (( ${#ASSEMBLIES[@]} == 0 )); then
    echo "No assemblies found in '$GENOMES_DIR' for the portals in '$PORTAL_LIST'" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUT_FILE")"
echo "Computing statistics for ${#ASSEMBLIES[@]} assemblies with $THREADS threads ..."

# Writes to a .part file first, so a failed run never overwrites an existing report.
# A 'portal' column (taken from the file name) is added in front of the seqkit columns.
# pipefail makes a seqkit error count as a failure, not only an awk error.
set -o pipefail
if ! seqkit stats --all --tabular --threads "$THREADS" "${ASSEMBLIES[@]}" \
    | awk -F '\t' -v OFS='\t' '
        NR == 1 { print "portal", $0; next }
        { portal = $1; sub(/.*\//, "", portal); sub(/_repeatmasked-assembly\.fasta\.gz$/, "", portal)
          print portal, $0 }' > "$OUT_FILE.part"; then
    echo "seqkit failed, report not written" >&2
    rm -f "$OUT_FILE.part"
    exit 1
fi
mv "$OUT_FILE.part" "$OUT_FILE"

echo "Saved statistics of ${#ASSEMBLIES[@]} assemblies to $OUT_FILE"
if (( ${#MISSING[@]} )); then
    echo "${#MISSING[@]} portals on the list have no assembly:" >&2
    printf '  %s\n' "${MISSING[@]}" >&2
fi
