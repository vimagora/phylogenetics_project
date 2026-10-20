#!/usr/bin/env python3
# select_busco_markers.py
# Description: Chooses the BUSCO markers for the species tree from the combined BUSCO full table
#              of the collapsed proteomes (combine_busco_tables.sh). A marker is kept when it is
#              complete and single-copy in most cohort genomes and rarely duplicated, since
#              duplicated markers may be paralogs. Run after clean_proteomes.py: every selected
#              sequence gets its clean ID, checked against the cohort ID map.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-20
# Version: 1.0
#
# Usage (from the project folder):
#   module load python-data
#   python3 scripts/select_busco_markers.py [--busco-table FILE] [--portal-list FILE] [--id-map FILE]
#       [--out-dir DIR] [--min-single-pct 90] [--max-dup-pct 5] [--max-genome-dup-pct 6]
#
# Rules, per marker, over the genomes on the portal list (the cohort; shadows never enter the tree):
#   - complete and single-copy in >= min-single-pct of the genomes;
#   - duplicated in <= max-dup-pct of the genomes, counting only genomes with <= max-genome-dup-pct
#     duplicated markers overall, so whole-genome duplications (e.g. Sacce1, Rhimi59_2) and
#     duplication-rich genomes (e.g. Neosp1) do not reject markers.
# In each genome, only a complete single-copy hit is used; duplicated and fragmented hits count
# as missing.
#
# Output, inside out-dir:
#   busco_marker_stats.tsv    counts per marker and whether it was selected
#   selected_markers.txt      selected BUSCO IDs, one per line
#   marker_sequences.tsv      selected marker, portal, BUSCO sequence name (first word of the collapsed
#                             proteome header) and clean ID (header in 05_clean_proteomes)
#   genome_coverage.tsv       selected markers present per genome, lowest first

import argparse
import os
import sys

import pandas as pd


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def read_portals(path):
    """Portal IDs from column 1 of a TSV with a header row."""
    with open(path) as f:
        lines = f.read().splitlines()[1:]
    portals = [line.split("\t")[0].strip() for line in lines]
    return [p for p in portals if p]


parser = argparse.ArgumentParser(description="Choose BUSCO markers for the species tree")
parser.add_argument("--busco-table", default="04_BUSCO/busco_full_table_collapsed_proteomes.tsv")
parser.add_argument("--portal-list", default="00_inputs/cohort_datasets.tsv")
parser.add_argument("--id-map", default="05_clean_proteomes/cohort/id_map.tsv")
parser.add_argument("--out-dir", default="06_species_tree")
parser.add_argument("--min-single-pct", type=float, default=90)
parser.add_argument("--max-dup-pct", type=float, default=5)
parser.add_argument("--max-genome-dup-pct", type=float, default=6)
args = parser.parse_args()

for path in (args.busco_table, args.portal_list, args.id_map):
    if not os.path.isfile(path):
        fail(f"'{path}' not found" + (", run scripts/clean_proteomes.py first" if path == args.id_map else ""))

portals = read_portals(args.portal_list)
if len(set(portals)) != len(portals):
    fail(f"duplicated portals in '{args.portal_list}'")

table = pd.read_csv(args.busco_table, sep="\t", dtype=str)
missing = sorted(set(portals) - set(table["portal"]))
if missing:
    fail(f"{len(missing)} portals have no BUSCO results in '{args.busco_table}': {', '.join(missing)}")
table = table[table["portal"].isin(portals)]

# One status per genome and marker (duplicated markers have one row per copy)
status = table.drop_duplicates(["portal", "Busco id"]).pivot(index="Busco id", columns="portal", values="Status")
if status.isna().any().any():
    fail("some genomes lack rows for some markers; were they run with different lineages?")

# Genomes with many duplicated markers do not count towards a marker's duplication rate
genome_dup_pct = (status == "Duplicated").mean() * 100
low_dup_genomes = genome_dup_pct.index[genome_dup_pct <= args.max_genome_dup_pct]
print(f"{len(portals)} genomes, {len(status)} markers; "
      f"{len(portals) - len(low_dup_genomes)} genomes with > {args.max_genome_dup_pct:g}% duplicated "
      f"markers are left out of the duplication rate: "
      + ", ".join(f"{p} ({genome_dup_pct[p]:.1f}%)" for p in genome_dup_pct.index.difference(low_dup_genomes)))

stats = pd.DataFrame({
    "n_single": (status == "Complete").sum(axis=1),
    "n_duplicated": (status == "Duplicated").sum(axis=1),
    "n_fragmented": (status == "Fragmented").sum(axis=1),
    "n_missing": (status == "Missing").sum(axis=1),
})
stats["pct_single"] = (stats["n_single"] / len(portals) * 100).round(1)
stats["pct_dup_low_dup_genomes"] = ((status[low_dup_genomes] == "Duplicated").mean(axis=1) * 100).round(1)
stats["selected"] = ((stats["pct_single"] >= args.min_single_pct)
                     & (stats["pct_dup_low_dup_genomes"] <= args.max_dup_pct))
descriptions = table.drop_duplicates("Busco id").set_index("Busco id")["Description"]
stats["description"] = descriptions.reindex(stats.index)
stats = stats.sort_values(["selected", "pct_single"], ascending=False)
selected = stats.index[stats["selected"]]
if len(selected) == 0:
    fail("no marker passes the thresholds")

# Protein ID of each selected marker in each genome where it is complete and single-copy
sequences = table[table["Busco id"].isin(selected) & (table["Status"] == "Complete")]
sequences = sequences[["Busco id", "portal", "Sequence"]].rename(columns={"Busco id": "busco_id", "Sequence": "busco_sequence"})
sequences = sequences.sort_values(["busco_id", "portal"])

# Clean ID from the BUSCO sequence name 'jgi|portal|proteinId|name', checked against the ID map
fields = sequences["busco_sequence"].str.split("|", expand=True).reindex(columns=range(3))
bad = (fields[0] != "jgi") | (fields[1] != sequences["portal"]) | ~fields[2].fillna("").str.fullmatch(r"\d+")
if bad.any():
    fail(f"{bad.sum()} BUSCO sequence names are not 'jgi|<portal>|<proteinId>|...', "
         f"e.g. {sequences.loc[bad, 'busco_sequence'].iloc[0]}")
sequences["clean_id"] = sequences["portal"] + "-" + fields[2]
clean_ids = set(pd.read_csv(args.id_map, sep="\t", usecols=["clean_id"], dtype=str)["clean_id"])
not_in_map = ~sequences["clean_id"].isin(clean_ids)
if not_in_map.any():
    fail(f"{not_in_map.sum()} selected sequences are not in '{args.id_map}', e.g. "
         f"{sequences.loc[not_in_map, 'clean_id'].iloc[0]}; were the clean proteomes made from the "
         f"same collapsed proteomes as the BUSCO run?")

coverage = sequences.groupby("portal").size().reindex(portals, fill_value=0).rename("n_markers").to_frame()
coverage["pct_markers"] = (coverage["n_markers"] / len(selected) * 100).round(1)
coverage["pct_duplicated_overall"] = genome_dup_pct.reindex(portals).round(1)
coverage = coverage.sort_values("n_markers")

os.makedirs(args.out_dir, exist_ok=True)
stats.to_csv(os.path.join(args.out_dir, "busco_marker_stats.tsv"), sep="\t", index_label="busco_id")
with open(os.path.join(args.out_dir, "selected_markers.txt"), "w") as f:
    f.write("\n".join(selected) + "\n")
sequences.to_csv(os.path.join(args.out_dir, "marker_sequences.tsv"), sep="\t", index=False)
coverage.to_csv(os.path.join(args.out_dir, "genome_coverage.tsv"), sep="\t", index_label="portal")

print(f"Selected {len(selected)} of {len(status)} markers (single-copy in >= {args.min_single_pct:g}% "
      f"of genomes, duplicated in <= {args.max_dup_pct:g}%), {len(sequences)} sequences")
print("Genomes with the fewest selected markers: "
      + ", ".join(f"{p} {row.pct_markers:g}%" for p, row in coverage.head(5).iterrows()))
print(f"Saved to {args.out_dir}/")
