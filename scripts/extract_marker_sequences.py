#!/usr/bin/env python3
# extract_marker_sequences.py
# Description: Writes one FASTA per selected BUSCO marker, with the marker's sequence from every
#              cohort genome where it is complete and single-copy (select_busco_markers.py).
#              Sequences come from the clean proteomes, and each is named by its portal only,
#              so the per-marker alignments can be concatenated and their trees compared.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-20
# Version: 1.0
#
# Usage (from the project folder):
#   module load python-data
#   python3 scripts/extract_marker_sequences.py [marker_sequences] [proteomes_dir] [output_dir]
#   marker_sequences  Table from select_busco_markers.py (default: 06_species_tree/marker_sequences.tsv).
#   proteomes_dir     Clean cohort proteomes (default: 05_clean_proteomes/cohort/proteomes).
#   output_dir        Where the FASTAs are saved (default: 06_species_tree/markers/raw).
#
# Output: <output_dir>/<busco_id>.faa, headers '>portal'. The clean ID of every sequence is in
# marker_sequences.tsv. Existing FASTAs are overwritten. The run stops if a clean ID is not found
# in its proteome, or if a marker has two sequences from one genome.

import csv
import os
import sys

table = sys.argv[1] if len(sys.argv) > 1 else "06_species_tree/marker_sequences.tsv"
proteomes_dir = sys.argv[2] if len(sys.argv) > 2 else "05_clean_proteomes/cohort/proteomes"
out_dir = sys.argv[3] if len(sys.argv) > 3 else "06_species_tree/markers/raw"
usage = f"Usage: {sys.argv[0]} [marker_sequences] [proteomes_dir] [output_dir]"

for path, kind in ((table, "file"), (proteomes_dir, "dir")):
    if not (os.path.isfile(path) if kind == "file" else os.path.isdir(path)):
        print(f"'{path}' not found\n{usage}", file=sys.stderr)
        sys.exit(1)

# clean_id -> (busco_id, portal), and the clean IDs wanted from each proteome
wanted, by_portal, markers = {}, {}, {}
with open(table) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        key = (row["busco_id"], row["portal"])
        if key in markers:
            print(f"ERROR: marker {key[0]} has two sequences from {key[1]}", file=sys.stderr)
            sys.exit(1)
        markers[key] = row["clean_id"]
        wanted[row["clean_id"]] = key
        by_portal.setdefault(row["portal"], set()).add(row["clean_id"])

# Read only the wanted sequences (clean proteomes have one sequence line per protein)
sequences = {}
for portal, ids in sorted(by_portal.items()):
    path = os.path.join(proteomes_dir, f"{portal}.faa")
    if not os.path.isfile(path):
        print(f"ERROR: proteome '{path}' not found", file=sys.stderr)
        sys.exit(1)
    with open(path) as f:
        clean_id = None
        for line in f:
            if line.startswith(">"):
                clean_id = line[1:].strip()
            elif clean_id in ids:
                sequences[clean_id] = line.strip()
    not_found = ids - sequences.keys()
    if not_found:
        print(f"ERROR ({portal}): {len(not_found)} clean IDs not in '{path}', e.g. {min(not_found)}",
              file=sys.stderr)
        sys.exit(1)

# One FASTA per marker, genomes in alphabetical order
by_marker = {}
for (busco_id, portal), clean_id in markers.items():
    by_marker.setdefault(busco_id, []).append((portal, sequences[clean_id]))
os.makedirs(out_dir, exist_ok=True)
for busco_id, records in by_marker.items():
    path = os.path.join(out_dir, f"{busco_id}.faa")
    with open(path + ".part", "w") as out:
        for portal, seq in sorted(records):
            out.write(f">{portal}\n{seq}\n")
    os.replace(path + ".part", path)

print(f"Wrote {len(by_marker)} marker FASTAs ({len(markers)} sequences from {len(by_portal)} genomes) "
      f"to {out_dir}")
