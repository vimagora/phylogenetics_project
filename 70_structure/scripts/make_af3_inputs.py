#!/usr/bin/env python3
# make_af3_inputs.py
# Family-specific: sugar transporters (structure concept test; see 70_structure/README.md).
# Description: Writes one AlphaFold3 input JSON per sugar transporter flagged for use in the family's
#              reference table (families/sugar_transporters/references.tsv), and a table linking each
#              structure ID (the reference ID) to its protein. Structure IDs are lowercase letters, digits and '_', because AlphaFold3
#              lowercases job names for its output folders.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-22
# Version: 1.1
#
# Usage (from the project folder):
#   python3 70_structure/scripts/make_af3_inputs.py [--out-dir DIR] [--single-sequence] [--only ID,ID,...]
#   --out-dir          Where the JSONs go (default: 70_structure/af3_inputs).
#   --single-sequence  Empty MSAs and no templates, so AlphaFold3 skips the database search. Only for
#                      quick tests of the container, weights and GPU; the models are much worse.
#   --only             Write only these structure IDs (e.g. for a one-protein test).
#
# Output: <out_dir>/<structure_id>.json, and 70_structure/structure_ids.tsv (all proteins, with
# source table, species, name, role, substrates, length, whether they are used, sequence, and the
# cohort portal and clean ID listed for the protein, where the reference table gives them).

import argparse
import csv
import json
import os
import sys

REFERENCES = "families/sugar_transporters/references.tsv"
IDS_TABLE = "70_structure/structure_ids.tsv"
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYX")


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


parser = argparse.ArgumentParser(description="Write AlphaFold3 input JSONs for the sugar transporters")
parser.add_argument("--out-dir", default="70_structure/af3_inputs")
parser.add_argument("--single-sequence", action="store_true")
parser.add_argument("--only", default="")
args = parser.parse_args()

if not os.path.isfile(REFERENCES):
    fail(f"'{REFERENCES}' not found")

# Structure IDs are the reference IDs; predicted proteins (the D. squalens transporters) are marked
# as such, and their substrates as predictions
proteins = []
with open(REFERENCES) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        predicted = row["evidence"] != "characterized"
        proteins.append(dict(id=row["reference_id"], source="dicsqu" if predicted else "characterized",
                             species=row["species"], name=row["protein_name"], role=row["role"],
                             substrates=f"predicted: {row['function']}" if predicted else row["function"],
                             sequence=row["sequence"], use=row["use"], portal=row["cohort_portal"],
                             clean_id=row["cohort_id"]))
for protein in proteins:
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
