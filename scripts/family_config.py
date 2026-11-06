#!/usr/bin/env python3
# family_config.py
# Description: Shared helpers of the generic family scripts: reads a family's settings
#              (families/<family>/family.conf) and the tables every family script needs.
#              Not run on its own; imported by summarize_family.py, build_family_fasta.py,
#              fetch_outgroups.py, make_family_itol.py, make_family_figure.py and build_subfamilies.py.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-11-05
# Version: 1.2

import csv
import os
import re
import sys

GENOMES = "00_inputs/genomes.tsv"
OG_DIR = "07_orthofinder/run/Results_cohort/Orthogroups"
PALETTE = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#1f78b4", "#b2df8a",
           "#fb9a99", "#cab2d6", "#ffff99", "#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69",
           "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f", "#9e0142", "#5e4fa2"]
MIN_COLOURED_ORTHOGROUP = 0.01  # share of a tree's tips an orthogroup needs for its own colour
# Phyla, largest in the cohort first, in an order checked for colour-blind separation; the largest avoid the
# blue of the coverage ring next to them. Further phyla are grey. Subphyla are shades of their phylum.
# Groups of non-fungal outgroups (e.g. plants) are greys.
TAXON_PALETTE = ["#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948", "#2a78d6", "#eb6834", "#1baf7a"]
EXTRA_PHYLUM = "#8c8c8c"
NO_SUBPHYLUM = "#e0e0e0"
NON_FUNGAL_GREYS = ["#525252", "#7f7f7f", "#a6a6a6", "#c8c8c8"]
UNKNOWN_TAXON = "#ececec"


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def read_config(family):
    """Settings of a family as a dict, plus its folder; stops if the family or a setting is missing."""
    folder = os.path.join("families", family)
    path = os.path.join(folder, "family.conf")
    if not os.path.isfile(path):
        fail(f"'{path}' not found; families are folders in families/")
    config = {"FAMILY": family, "FAMILY_DIR": folder}
    with open(path) as f:
        for line in f:
            match = re.match(r"\s*([A-Z_]+)=(.*)$", line)
            if match and not line.lstrip().startswith("#"):
                config[match.group(1)] = match.group(2).strip().strip('"')
    for key in ("NAME", "HMMS", "MIN_HMM_COVERAGE", "OUT_DIR"):
        if not config.get(key):
            fail(f"{key} is not set in {path}")
    config["HMMS"] = config["HMMS"].split()
    config["MIN_HMM_COVERAGE"] = float(config["MIN_HMM_COVERAGE"])
    return config


def read_tsv(path, required=True):
    if not os.path.isfile(path):
        if required:
            fail(f"'{path}' not found")
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, header, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def taxonomy():
    """Taxonomy of the cohort genomes (role = cohort in genomes.tsv) by portal; shadows are left out."""
    return {r["portal"].strip(): {k: (v or "").strip() for k, v in r.items()} for r in read_tsv(GENOMES)
            if r["role"].strip() == "cohort"}


def lineage_of(family_dir):
    """A function giving a tip's phylum, subphylum or class (rank): from genomes.tsv for cohort proteins,
    orphans.tsv for orphans, and the outgroup's group (e.g. plants) for outgroups; 'unknown' otherwise.
    tip: a row of a tips.tsv."""
    tax = taxonomy()
    orphans = {r["pseudo_portal"]: r for r in read_tsv(os.path.join(family_dir, "orphans.tsv"), required=False)}
    outgroups = {r["outgroup_id"]: r for r in read_tsv(os.path.join(family_dir, "outgroups.tsv"), required=False)}

    def lineage(tip, rank):
        if tip["kind"] == "outgroup" and tip["tip_id"] in outgroups:
            return outgroups[tip["tip_id"]]["group"]
        row = tax.get(tip["portal"]) or orphans.get(tip["portal"])
        return (row[rank] or "unknown") if row else "unknown"
    return lineage


def shade(colour, t):
    """A colour darkened (t < 0) or lightened (t > 0) by the share t."""
    rgb = [int(colour[i:i + 2], 16) for i in (1, 3, 5)]
    target = 0 if t < 0 else 255
    return "#" + "".join(f"{round(c + (target - c) * abs(t)):02x}" for c in rgb)


def taxon_colours(family_dir):
    """Colours of phyla and subphyla, the same in every figure and iTOL file of a family: phyla in the order
    of their genome count in the cohort, subphyla as shades of their phylum (darkest for the largest), and
    the groups of non-fungal outgroups in greys. Returns (phylum colours, subphylum colours)."""
    counts, subphyla = {}, {}
    rows = list(taxonomy().values()) + read_tsv(os.path.join(family_dir, "orphans.tsv"), required=False)
    for row in rows:
        counts[row["phylum"]] = counts.get(row["phylum"], 0) + 1
        if row["subphylum"] and row["subphylum"] != "none":
            sub = subphyla.setdefault(row["phylum"], {})
            sub[row["subphylum"]] = sub.get(row["subphylum"], 0) + 1
    phyla = sorted((p for p in counts if p), key=lambda p: (-counts[p], p))
    phylum_colours = {p: TAXON_PALETTE[i] if i < len(TAXON_PALETTE) else EXTRA_PHYLUM for i, p in enumerate(phyla)}
    subphylum_colours = {}
    for phylum, sizes in subphyla.items():
        names = sorted(sizes, key=lambda s: (-sizes[s], s))
        for i, name in enumerate(names):
            t = -0.3 + 0.8 * i / (len(names) - 1) if len(names) > 1 else 0.0
            subphylum_colours[name] = shade(phylum_colours[phylum], t)
    groups = sorted({r["group"] for r in read_tsv(os.path.join(family_dir, "outgroups.tsv"), required=False)})
    for i, group in enumerate(groups):
        phylum_colours[group] = subphylum_colours[group] = NON_FUNGAL_GREYS[min(i, len(NON_FUNGAL_GREYS) - 1)]
    return phylum_colours, subphylum_colours


def portal_of(protein_id):
    """Portal (or pseudo-portal) of a '<portal>-<protein>' ID: everything before the last '-'."""
    return protein_id.rsplit("-", 1)[0]


def read_fasta(path):
    sequences, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = []
            elif name:
                sequences[name].append(line)
    return {k: "".join(v) for k, v in sequences.items()}


def read_domtbl(path):
    """Hits of an hmmsearch domain table by protein: full-sequence score, and coverage of the HMM as the
    share of its positions covered by the protein's domains, so that the pieces of a split domain add up
    and repeated domains count once."""
    proteins = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.split()
            p = proteins.setdefault(x[0], dict(length=int(x[2]), hmm=x[4], hmm_length=int(x[5]),
                                               evalue=float(x[6]), score=float(x[7]), spans=[]))
            p["spans"].append((int(x[15]), int(x[16])))
    for p in proteins.values():
        covered = set()
        for a, b in p["spans"]:
            covered.update(range(a, b + 1))
        p["coverage"] = len(covered) / p["hmm_length"]
    return proteins


def orthogroup_colours(tips):
    """Colours of the orthogroups holding at least 1% of a tree's tips, largest first, up to the size of
    the palette; the iTOL files and the figure share them. tips: rows of tree/tips.tsv."""
    sizes = {}
    for t in tips:
        if t["orthogroup"]:
            sizes[t["orthogroup"]] = sizes.get(t["orthogroup"], 0) + 1
    large = [og for og in sorted(sizes, key=lambda og: (-sizes[og], og))
             if sizes[og] >= MIN_COLOURED_ORTHOGROUP * len(tips)]
    return dict(zip(large, PALETTE))
