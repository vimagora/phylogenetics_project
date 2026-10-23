#!/usr/bin/env python3
# make_family_itol.py
# Description: Writes iTOL annotation files for a gene family's tree (build_family_fasta.py tips),
#              for review by eye. Function is shown only for characterized references (not for
#              predictions, sensors or unanchored members); nothing is inferred for other tips.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-23
# Version: 1.0
#
# Usage (from the project folder; standard library only): python3 scripts/make_family_itol.py <family>
#
# Reads <OUT_DIR>/tree/tips.tsv, <OUT_DIR>/anchors.tsv, and from families/<family>/: references.tsv,
# orphans.tsv, outgroups.tsv and, if present, function_columns.tsv (label, colour, terms: the words in
# a reference's function that tick that column, comma-separated, lowercase).
#
# Output, inside <OUT_DIR>/itol (drag the files onto the tree in iTOL):
#   labels.txt      tip labels: ID, reference names and species
#   kind.txt        colour strip: characterized anchor, probable anchor, predicted, sensor, orphan,
#                   outgroup, member
#   phylum.txt, class.txt   colour strips of the tip's genome (orphans: orphans.tsv; outgroups: group)
#   orthogroup.txt  colour strip of the 12 largest orthogroups (others grey)
#   function.txt    binary dataset of function_columns.tsv, characterized references only
#                   (without function_columns.tsv: function as text, characterized references only)

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_config import fail, read_config, read_tsv, taxonomy  # noqa: E402

PALETTE = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#1f78b4", "#b2df8a",
           "#fb9a99", "#cab2d6", "#ffff99", "#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69",
           "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f", "#9e0142", "#5e4fa2"]
KIND_COLOURS = {"characterized anchor": "#e41a1c", "probable anchor": "#ff7f00", "predicted (D. squalens)": "#377eb8",
                "sensor": "#984ea3", "orphan (characterized)": "#a65628", "outgroup": "#000000", "member": "#dddddd"}

if len(sys.argv) != 2:
    fail("usage: python3 scripts/make_family_itol.py <family>")
config = read_config(sys.argv[1])
out = config["OUT_DIR"]
tips = read_tsv(os.path.join(out, "tree", "tips.tsv"))
references = {r["reference_id"]: r for r in read_tsv(os.path.join(config["FAMILY_DIR"], "references.tsv"))}
orphans = {r["pseudo_portal"]: r for r in read_tsv(os.path.join(config["FAMILY_DIR"], "orphans.tsv"))}
outgroups = {r["outgroup_id"]: r for r in read_tsv(os.path.join(config["FAMILY_DIR"], "outgroups.tsv"))}
columns = read_tsv(os.path.join(config["FAMILY_DIR"], "function_columns.tsv"), required=False)
tax = taxonomy()


def refs_of(tip):
    return [references[r] for r in tip["references"].split(",") if r]


def kind(tip):
    if tip["kind"] == "outgroup":
        return "outgroup"
    refs = refs_of(tip)
    if not refs:
        return "member"
    if any(r["role"] == "sensor" for r in refs):
        return "sensor"
    if tip["kind"] == "orphan":
        return "orphan (characterized)" if refs[0]["evidence"] == "characterized" else "predicted (D. squalens)"
    if all(r["evidence"] != "characterized" for r in refs):
        return "predicted (D. squalens)"
    return "probable anchor" if "probable_anchor" in tip["status"] else "characterized anchor"


def lineage(tip, rank):
    if tip["kind"] == "outgroup":
        return outgroups[tip["tip_id"]]["group"]
    row = tax.get(tip["portal"]) or orphans.get(tip["portal"])
    return row[rank] if row else "unknown"


def functions(tip):
    """Function words of the tip's characterized, transporting references."""
    words = set()
    for r in refs_of(tip):
        if r["evidence"] == "characterized" and r["role"] != "sensor":
            words |= {w.strip().lower() for w in re.split(r"[;,]", r["function"]) if w.strip()}
    return words


def species(tip):
    if tip["kind"] == "outgroup":
        return outgroups[tip["tip_id"]]["organism"]
    row = tax.get(tip["portal"]) or orphans.get(tip["portal"])
    return " ".join(row["species"].split()[:2]) if row else tip["portal"]


def colour_strip(name, label, values, legend=True):
    groups = sorted(set(values.values()), key=lambda g: (g in ("unknown", "other"), g))
    colours = {g: ("#dddddd" if g in ("unknown", "other") else PALETTE[i % len(PALETTE)]) for i, g in enumerate(groups)}
    if name == "kind.txt":
        colours = KIND_COLOURS
        groups = [g for g in KIND_COLOURS if g in set(values.values())]
    lines = ["DATASET_COLORSTRIP", "SEPARATOR TAB", f"DATASET_LABEL\t{label}", "COLOR\t#000000"]
    if legend:
        lines += [f"LEGEND_TITLE\t{label}", "LEGEND_SHAPES\t" + "\t".join("1" for _ in groups),
                  "LEGEND_COLORS\t" + "\t".join(colours[g] for g in groups), "LEGEND_LABELS\t" + "\t".join(groups)]
    lines.append("DATA")
    lines += [f"{t}\t{colours[v]}\t{v}" for t, v in values.items()]
    with open(os.path.join(itol, name), "w") as f:
        f.write("\n".join(lines) + "\n")


itol = os.path.join(out, "itol")
os.makedirs(itol, exist_ok=True)

with open(os.path.join(itol, "labels.txt"), "w") as f:
    f.write("LABELS\nSEPARATOR TAB\nDATA\n")
    for t in tips:
        names = "/".join(r["protein_name"] for r in refs_of(t))
        if t["kind"] == "outgroup":
            names = outgroups[t["tip_id"]]["protein"]
        f.write(f"{t['tip_id']}\t{t['tip_id']}{' ' + names if names else ''} ({species(t)})\n")

colour_strip("kind.txt", "Kind", {t["tip_id"]: kind(t) for t in tips})
colour_strip("phylum.txt", "Phylum", {t["tip_id"]: lineage(t, "phylum") for t in tips})
colour_strip("class.txt", "Class", {t["tip_id"]: lineage(t, "class") for t in tips}, legend=False)
sizes = {}
for t in tips:
    if t["orthogroup"]:
        sizes[t["orthogroup"]] = sizes.get(t["orthogroup"], 0) + 1
top = set(sorted(sizes, key=lambda og: -sizes[og])[:12])
colour_strip("orthogroup.txt", "Orthogroup",
             {t["tip_id"]: (t["orthogroup"] if t["orthogroup"] in top else "other") for t in tips})

characterized = {t["tip_id"]: functions(t) for t in tips if functions(t)}
if columns:
    lines = ["DATASET_BINARY", "SEPARATOR TAB", "DATASET_LABEL\tFunction (characterized)", "COLOR\t#000000",
             "FIELD_SHAPES\t" + "\t".join("1" for _ in columns), "FIELD_LABELS\t" + "\t".join(c["label"] for c in columns),
             "FIELD_COLORS\t" + "\t".join(c["colour"] for c in columns), "SHOW_LABELS\t1", "DATA"]
    terms = [{w.strip() for w in c["terms"].split(",")} for c in columns]
    for tip, words in characterized.items():
        lines.append(tip + "\t" + "\t".join("1" if words & ts else "-1" for ts in terms))
else:
    lines = ["DATASET_TEXT", "SEPARATOR TAB", "DATASET_LABEL\tFunction (characterized)", "COLOR\t#000000", "DATA"]
    lines += [f"{tip}\t{', '.join(sorted(words))}\t-1\t#000000\tnormal\t1" for tip, words in characterized.items()]
with open(os.path.join(itol, "function.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")

counts = {}
for t in tips:
    counts[kind(t)] = counts.get(kind(t), 0) + 1
unmatched = sorted({w for words in characterized.values() for w in words}
                   - set().union(*[{x.strip() for x in c["terms"].split(",")} for c in columns])) if columns else []
print(f"iTOL files for {len(tips)} tips in {itol}/: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
if unmatched:
    print(f"  function words in no column of function_columns.tsv: {', '.join(unmatched)}")
