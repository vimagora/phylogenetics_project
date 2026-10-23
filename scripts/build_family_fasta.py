#!/usr/bin/env python3
# build_family_fasta.py
# Description: Collects the sequences of a gene family's tree (summarize_family.py results):
#              - members: cohort proteins whose family HMM coverage passes MIN_HMM_COVERAGE, as full
#                sequences from the clean proteomes (anchored references are among them);
#              - orphans: references without a cohort copy, as '<pseudo-portal>-<protein>', with the
#                reference's sequence, if the reference itself passes the membership rule;
#              - outgroups: fetch_outgroups.py.
#              Shadows never enter. Run where the clean proteomes are (Roihu).
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-23
# Version: 1.0
#
# Usage (from the project folder; standard library only):
#   python3 scripts/build_family_fasta.py <family> [proteomes_dir]
#   proteomes_dir  Clean cohort proteomes (default: 05_clean_proteomes/cohort/proteomes).
#
# Output, inside <OUT_DIR>/tree:
#   sequences.faa   one sequence per tip, unaligned
#   tips.tsv        tip_id, kind (member, orphan, outgroup), portal (pseudo-portal for orphans,
#                   abbreviation for outgroups), orthogroup, references (anchored or orphan reference
#                   IDs) and anchor status

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_config import fail, portal_of, read_config, read_fasta, read_tsv, write_tsv  # noqa: E402

if len(sys.argv) not in (2, 3):
    fail("usage: python3 scripts/build_family_fasta.py <family> [proteomes_dir]")
config = read_config(sys.argv[1])
proteomes_dir = sys.argv[2] if len(sys.argv) == 3 else "05_clean_proteomes/cohort/proteomes"
out = config["OUT_DIR"]
members = read_tsv(os.path.join(out, "members.tsv"))
anchors = read_tsv(os.path.join(out, "anchors.tsv"))
references = {r["reference_id"]: r for r in read_tsv(os.path.join(config["FAMILY_DIR"], "references.tsv"))}
orphan_table = {r["pseudo_portal"] for r in read_tsv(os.path.join(config["FAMILY_DIR"], "orphans.tsv"))}
outgroups_faa = os.path.join(out, "outgroups.faa")
if not os.path.isfile(outgroups_faa):
    fail(f"'{outgroups_faa}' not found; run scripts/fetch_outgroups.py {config['FAMILY']}")

sequences, tips = {}, []

# Members, from the clean proteomes of their genomes
by_portal = {}
for m in members:
    by_portal.setdefault(m["portal"], []).append(m)
for portal, rows in sorted(by_portal.items()):
    path = os.path.join(proteomes_dir, f"{portal}.faa")
    if not os.path.isfile(path):
        fail(f"'{path}' not found")
    proteome = read_fasta(path)
    for m in rows:
        if m["clean_id"] not in proteome:
            fail(f"{m['clean_id']} is not in {path}")
        sequences[m["clean_id"]] = proteome[m["clean_id"]].rstrip("*")
        tips.append([m["clean_id"], "member", portal, m["orthogroup"], m["anchor_of"], m["anchor_status"]])

# Orphans that are family members themselves
skipped = []
for a in anchors:
    if not a["status"].startswith("orphan"):
        continue
    if float(a["reference_hmm_coverage"] or 0) < config["MIN_HMM_COVERAGE"]:
        skipped.append(a["reference_id"])
        continue
    pseudo = portal_of(a["tree_id"])
    if pseudo not in orphan_table:
        fail(f"pseudo-portal {pseudo} ({a['reference_id']}) is not described in families/{config['FAMILY']}/orphans.tsv")
    sequences[a["tree_id"]] = references[a["reference_id"]]["sequence"]
    tips.append([a["tree_id"], "orphan", pseudo, "", a["reference_id"], a["status"]])

# Outgroups
for name, seq in read_fasta(outgroups_faa).items():
    sequences[name] = seq
    tips.append([name, "outgroup", portal_of(name), "", "", ""])

if len(sequences) != len(tips):
    fail("two tips have the same ID")
tree_dir = os.path.join(out, "tree")
os.makedirs(tree_dir, exist_ok=True)
with open(os.path.join(tree_dir, "sequences.faa"), "w") as f:
    for tip in tips:
        f.write(f">{tip[0]}\n{sequences[tip[0]]}\n")
write_tsv(os.path.join(tree_dir, "tips.tsv"), ["tip_id", "kind", "portal", "orthogroup", "references", "status"], tips)
kinds = {k: sum(t[1] == k for t in tips) for k in ("member", "orphan", "outgroup")}
print(f"{len(tips)} sequences: {kinds['member']} members, {kinds['orphan']} orphans, {kinds['outgroup']} outgroups"
      + (f"; left out (below the membership rule): {', '.join(skipped)}" if skipped else "")
      + f". Saved to {tree_dir}/")
