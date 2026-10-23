#!/usr/bin/env python3
# fetch_outgroups.py
# Description: Downloads the outgroup sequences of a gene family from UniProt, as listed in
#              families/<family>/outgroups.tsv (outgroup_id, group, organism, protein, uniprot,
#              root, notes). Outgroups form their own branch, whose expected internal structure
#              (families/<family>/expectations.md) is a check on alignment and tree building;
#              those with root = TRUE are used to root the tree.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-23
# Version: 1.0
#
# Usage (from the project folder; needs internet): python3 scripts/fetch_outgroups.py <family>
#
# Output: <OUT_DIR>/outgroups.faa, headers '>outgroup_id' (OUT_DIR from family.conf).

import os
import re
import sys
import time
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_config import fail, read_config, read_tsv  # noqa: E402

if len(sys.argv) != 2:
    fail("usage: python3 scripts/fetch_outgroups.py <family>")
config = read_config(sys.argv[1])
outgroups = read_tsv(os.path.join(config["FAMILY_DIR"], "outgroups.tsv"))
ids = [o["outgroup_id"] for o in outgroups]
if len(ids) != len(set(ids)):
    fail("outgroup IDs in outgroups.tsv are not unique")
bad = [i for i in ids if not re.fullmatch(r"[A-Za-z0-9_]+-[A-Za-z0-9_.]+", i)]
if bad:
    fail(f"outgroup IDs must look like '<Abbrev>-<protein>' (letters, digits, '_'): {', '.join(bad)}")

records = []
for o in outgroups:
    for attempt in range(3):
        try:
            with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{o['uniprot']}.fasta", timeout=60) as r:
                text = r.read().decode()
            break
        except OSError:
            time.sleep(5 * (attempt + 1))
    else:
        fail(f"could not download {o['uniprot']}")
    header, *lines = text.strip().splitlines()
    if f"|{o['uniprot']}|" not in header:
        fail(f"UniProt returned {header[:60]} for {o['uniprot']}")
    records.append((o["outgroup_id"], "".join(lines)))
    print(f"  {o['outgroup_id']}: {o['uniprot']}, {len(records[-1][1])} aa")

os.makedirs(config["OUT_DIR"], exist_ok=True)
out = os.path.join(config["OUT_DIR"], "outgroups.faa")
with open(out, "w") as f:
    for name, seq in records:
        f.write(f">{name}\n{seq}\n")
print(f"{len(records)} outgroups ({sum(o['root'] == 'TRUE' for o in outgroups)} for rooting) saved to {out}")
