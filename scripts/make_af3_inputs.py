#!/usr/bin/env python3
# make_af3_inputs.py
# Family-specific: sugar transporters (structure concept test; see README, "Sugar transporter structures").
# Description: Writes one AlphaFold3 input JSON per sugar transporter flagged for structure analysis
#              in the characterized and D. squalens tables, and a table linking each structure ID to
#              its protein. Structure IDs are lowercase letters, digits and '_', because AlphaFold3
#              lowercases job names for its output folders.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-22
# Version: 1.0
#
# Usage (from the project folder):
#   python3 scripts/make_af3_inputs.py [--out-dir DIR] [--single-sequence] [--only ID,ID,...]
#   --out-dir          Where the JSONs go (default: 08_structures/af3_inputs).
#   --single-sequence  Empty MSAs and no templates, so AlphaFold3 skips the database search. Only for
#                      quick tests of the container, weights and GPU; the models are much worse.
#   --only             Write only these structure IDs (e.g. for a one-protein test).
#
# Output: <out_dir>/<structure_id>.json, and 08_structures/structure_ids.tsv (all proteins, with
# source table, species, name, role, substrates, length, whether they are used, sequence, and the
# cohort portal and clean ID listed for the protein, where the input tables give them).

import argparse
import csv
import json
import os
import re
import sys

CHARACTERIZED = "00_inputs/characterized_sugar_transporters.tsv"
DICSQU = "00_inputs/dicsqu_sugar_transporters.tsv"
IDS_TABLE = "08_structures/structure_ids.tsv"
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYX")


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def slug(text):
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", text.lower())).strip("_")


parser = argparse.ArgumentParser(description="Write AlphaFold3 input JSONs for the sugar transporters")
parser.add_argument("--out-dir", default="08_structures/af3_inputs")
parser.add_argument("--single-sequence", action="store_true")
parser.add_argument("--only", default="")
args = parser.parse_args()

for path in (CHARACTERIZED, DICSQU):
    if not os.path.isfile(path):
        fail(f"'{path}' not found")

proteins = []
with open(CHARACTERIZED) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        listed = re.search(r"clean ID (\S+?)[;,)]?(\s|$)", row["notes"])
        proteins.append(dict(base=slug(f"{row['short_name']}_{row['protein_name']}"), source="characterized",
                             species=row["species"], name=row["protein_name"], role=row["role"],
                             substrates=row["group_a"], sequence=row["sequence"],
                             use=row["use_in_structure_analysis"], portal=row["portal_id"].strip(),
                             clean_id=listed.group(1) if listed else ""))
with open(DICSQU) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        proteins.append(dict(base=slug(f"dicsqu_{row['gene']}"), source="dicsqu",
                             species="Dichomitus squalens", name=row["gene"], role="transporter",
                             substrates=f"predicted: {row['predicted_specificity']}",
                             sequence=row["sequence"], use=row["use_in_structure_analysis"],
                             portal=row["clean_id"].rsplit("-", 1)[0], clean_id=row["clean_id"]))

# Unique structure IDs: a numbered suffix only where two proteins would get the same ID
counts = {}
for protein in proteins:
    counts[protein["base"]] = counts.get(protein["base"], 0) + 1
seen = {}
for protein in proteins:
    base = protein["base"]
    seen[base] = seen.get(base, 0) + 1
    protein["id"] = base if counts[base] == 1 else f"{base}_{seen[base]}"
    odd = set(protein["sequence"]) - AMINO_ACIDS
    if odd:
        fail(f"{protein['id']} has characters other than amino acids: {''.join(sorted(odd))}")

only = {x.strip() for x in args.only.split(",") if x.strip()}
unknown = only - {p["id"] for p in proteins}
if unknown:
    fail(f"unknown structure IDs in --only: {', '.join(sorted(unknown))}")

os.makedirs(args.out_dir, exist_ok=True)
os.makedirs(os.path.dirname(IDS_TABLE), exist_ok=True)
n_written = 0
for protein in proteins:
    if protein["use"] != "TRUE" or (only and protein["id"] not in only):
        continue
    chain = {"id": "A", "sequence": protein["sequence"]}
    if args.single_sequence:
        chain.update(unpairedMsa="", pairedMsa="", templates=[])
    job = {"name": protein["id"], "modelSeeds": [1], "sequences": [{"protein": chain}],
           "dialect": "alphafold3", "version": 1}
    with open(os.path.join(args.out_dir, f"{protein['id']}.json"), "w") as f:
        json.dump(job, f, indent=1)
    n_written += 1

with open(IDS_TABLE, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
    writer.writerow(["structure_id", "source", "species", "protein_name", "role", "substrates",
                     "length", "use_in_structure_analysis", "sequence", "cohort_portal", "cohort_id"])
    for p in proteins:
        writer.writerow([p["id"], p["source"], p["species"], p["name"], p["role"], p["substrates"],
                         len(p["sequence"]), p["use"], p["sequence"], p["portal"], p["clean_id"]])

mode = "single-sequence (test only)" if args.single_sequence else "with MSA search"
print(f"Wrote {n_written} AlphaFold3 inputs ({mode}) to {args.out_dir}; IDs in {IDS_TABLE}")
