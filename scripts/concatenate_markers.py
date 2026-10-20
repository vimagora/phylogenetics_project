#!/usr/bin/env python3
# concatenate_markers.py
# Description: Filters the trimmed marker alignments (align_markers.sbatch) and concatenates the
#              kept ones into a supermatrix with a partition file, one partition per marker, for
#              the concatenated IQ-TREE run. The kept markers are also the ones used for the
#              per-marker trees.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-20
# Version: 1.0
#
# Usage (from the project folder):
#   module load python-data
#   python3 scripts/concatenate_markers.py [--trimmed-dir DIR] [--portal-list FILE] [--out-dir DIR]
#       [--min-length 100] [--min-taxa-pct 90]
#
# A sequence that is only gaps after trimming counts as missing. A marker is kept when its trimmed
# alignment has >= min-length columns and sequences from >= min-taxa-pct of the cohort genomes.
# Genomes missing from a marker are filled with gaps in the supermatrix.
#
# Output, inside out-dir:
#   supermatrix.faa        concatenated alignment, one sequence per cohort genome (portal names)
#   partitions.nex         NEXUS charset per kept marker, for 'iqtree3 -p'
#   kept_markers.txt       kept BUSCO IDs, one per line
#   marker_summary.tsv     length, genomes and gap % of every trimmed marker, and whether it was kept
#   taxon_occupancy.tsv    kept markers and % of supermatrix columns with data, per genome, lowest first

import argparse
import os
import sys


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def read_portals(path):
    """Portal IDs from column 1 of a TSV with a header row."""
    with open(path) as f:
        lines = f.read().splitlines()[1:]
    portals = [line.split("\t")[0].strip() for line in lines]
    return [p for p in portals if p]


def read_alignment(path):
    """{name: aligned sequence}; stops if names repeat or lengths differ."""
    records, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in records:
                    fail(f"'{path}' has two sequences named {name}")
                records[name] = []
            elif line:
                records[name].append(line)
    records = {n: "".join(s).upper() for n, s in records.items()}
    if len({len(s) for s in records.values()}) > 1:
        fail(f"'{path}' is not an alignment (sequences of different lengths)")
    return records


parser = argparse.ArgumentParser(description="Filter and concatenate trimmed marker alignments")
parser.add_argument("--trimmed-dir", default="06_species_tree/markers/trimmed")
parser.add_argument("--portal-list", default="00_inputs/cohort_datasets.tsv")
parser.add_argument("--out-dir", default="06_species_tree/concatenated")
parser.add_argument("--min-length", type=int, default=100)
parser.add_argument("--min-taxa-pct", type=float, default=90)
args = parser.parse_args()

for path in (args.trimmed_dir, args.portal_list):
    if not os.path.exists(path):
        fail(f"'{path}' not found")

portals = read_portals(args.portal_list)
files = sorted(f for f in os.listdir(args.trimmed_dir) if f.endswith(".faa"))
if not files:
    fail(f"no trimmed alignments in '{args.trimmed_dir}', run scripts/align_markers.sbatch first")

summary, kept = [], []
for file in files:
    busco_id = file.removesuffix(".faa")
    records = read_alignment(os.path.join(args.trimmed_dir, file))
    unknown = records.keys() - set(portals)
    if unknown:
        fail(f"'{file}' has sequences not on the portal list: {', '.join(sorted(unknown))}")
    records = {n: s for n, s in records.items() if s.strip("-X?")}
    length = len(next(iter(records.values()), ""))
    n_taxa = len(records)
    gap_pct = (100 * sum(s.count("-") for s in records.values()) / (length * n_taxa)) if n_taxa and length else 100
    keep = length >= args.min_length and n_taxa >= args.min_taxa_pct / 100 * len(portals)
    summary.append((busco_id, length, n_taxa, round(gap_pct, 1), keep))
    if keep:
        kept.append((busco_id, length, records))

if not kept:
    fail("no marker passes the filters")

os.makedirs(args.out_dir, exist_ok=True)
with open(os.path.join(args.out_dir, "marker_summary.tsv"), "w") as f:
    f.write("busco_id\ttrimmed_length\tn_genomes\tpct_gaps\tkept\n")
    for row in summary:
        f.write("\t".join(map(str, row)) + "\n")
with open(os.path.join(args.out_dir, "kept_markers.txt"), "w") as f:
    f.write("".join(f"{busco_id}\n" for busco_id, _, _ in kept))

# Supermatrix and partitions; genomes missing from a marker get gaps
total = sum(length for _, length, _ in kept)
with open(os.path.join(args.out_dir, "partitions.nex"), "w") as f:
    f.write("#nexus\nbegin sets;\n")
    start = 1
    for busco_id, length, _ in kept:
        f.write(f"    charset {busco_id} = {start}-{start + length - 1};\n")
        start += length
    f.write("end;\n")
occupancy = []
with open(os.path.join(args.out_dir, "supermatrix.faa.part"), "w") as f:
    for portal in portals:
        parts = [records.get(portal, "-" * length) for _, length, records in kept]
        seq = "".join(parts)
        n_markers = sum(portal in records for _, _, records in kept)
        data_pct = 100 * sum(c not in "-X?" for c in seq) / total
        occupancy.append((portal, n_markers, round(data_pct, 1)))
        f.write(f">{portal}\n{seq}\n")
os.replace(os.path.join(args.out_dir, "supermatrix.faa.part"), os.path.join(args.out_dir, "supermatrix.faa"))
occupancy.sort(key=lambda row: row[2])
with open(os.path.join(args.out_dir, "taxon_occupancy.tsv"), "w") as f:
    f.write("portal\tn_markers\tpct_columns_with_data\n")
    for row in occupancy:
        f.write("\t".join(map(str, row)) + "\n")

print(f"Kept {len(kept)} of {len(files)} markers (>= {args.min_length} columns after trimming, "
      f">= {args.min_taxa_pct:g}% of {len(portals)} genomes): supermatrix of {total} columns")
print("Genomes with the least data: "
      + ", ".join(f"{p} {pct:g}% ({n} markers)" for p, n, pct in occupancy[:5]))
print(f"Saved to {args.out_dir}/")
