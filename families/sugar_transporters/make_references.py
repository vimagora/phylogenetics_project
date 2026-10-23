#!/usr/bin/env python3
# make_references.py
# Description: Family-specific (sugar transporters). One-off converter of the two curated sugar
#              transporter tables in 00_inputs into the generic reference table of the family
#              scripts, families/sugar_transporters/references.tsv. Reference IDs are the same as
#              the structure IDs of the concept test (make_af3_inputs.py). New families write their
#              reference table directly in the generic format (see README, "Gene families").
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-23
# Version: 1.0
#
# Usage (from the project folder): python3 families/sugar_transporters/make_references.py
#
# Output: families/sugar_transporters/references.tsv, one row per protein, with
#   reference_id   lowercase ID (letters, digits, '_')
#   evidence       characterized (literature) or predicted (D. squalens predictions)
#   species, short_name, protein_name, role (transporter or sensor)
#   function       substrates (';'-separated for characterized proteins)
#   cohort_portal  the protein's own cohort genome, where the curated table gives it
#   cohort_id      its clean ID in that genome, where known
#   use            TRUE if it enters the family analysis
#   references, notes, sequence

import csv
import os
import re
import sys

CHARACTERIZED = "00_inputs/characterized_sugar_transporters.tsv"
DICSQU = "00_inputs/dicsqu_sugar_transporters.tsv"
OUT = "families/sugar_transporters/references.tsv"
COLUMNS = ["reference_id", "evidence", "species", "short_name", "protein_name", "role", "function", "cohort_portal",
           "cohort_id", "use", "references", "notes", "sequence"]


def slug(text):
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", text.lower())).strip("_")


for path in (CHARACTERIZED, DICSQU):
    if not os.path.isfile(path):
        sys.exit(f"ERROR: '{path}' not found")

rows = []
with open(CHARACTERIZED) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        listed = re.search(r"clean ID (\S+?)[;,)]?(\s|$)", r["notes"])
        rows.append(dict(base=slug(f"{r['short_name']}_{r['protein_name']}"), evidence="characterized",
                         species=r["species"], short_name=r["short_name"], protein_name=r["protein_name"],
                         role=r["role"], function=r["group_a"], cohort_portal=r["portal_id"].strip(),
                         cohort_id=listed.group(1) if listed else "", use=r["use_in_structure_analysis"],
                         references=r["references"], notes=r["notes"], sequence=r["sequence"]))
with open(DICSQU) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows.append(dict(base=slug(f"dicsqu_{r['gene']}"), evidence="predicted", species="Dichomitus squalens",
                         short_name="Dicsqu", protein_name=r["gene"], role="transporter",
                         function=r["predicted_specificity"], cohort_portal=r["clean_id"].rsplit("-", 1)[0],
                         cohort_id=r["clean_id"], use=r["use_in_structure_analysis"], references="",
                         notes=r["notes"], sequence=r["sequence"]))

# Unique IDs: a numbered suffix only where two proteins would get the same ID (as in make_af3_inputs.py)
counts, seen = {}, {}
for r in rows:
    counts[r["base"]] = counts.get(r["base"], 0) + 1
for r in rows:
    seen[r["base"]] = seen.get(r["base"], 0) + 1
    r["reference_id"] = r["base"] if counts[r["base"]] == 1 else f"{r['base']}_{seen[r['base']]}"

with open(OUT, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
print(f"{len(rows)} references ({sum(r['use'] == 'TRUE' for r in rows)} used) written to {OUT}")
