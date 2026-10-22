#!/usr/bin/env python3
# make_itol_taxonomy.py
# Description: Writes iTOL colour-strip annotation files (one per rank: phylum, subphylum, class,
#              order) for a species tree, from the cohort taxonomy table. Drag the files onto the
#              tree in iTOL. Colours follow the order of the tips in the tree, so neighbouring
#              groups get different colours; each group keeps one colour even if it is split.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-21
# Version: 1.0
#
# Usage (from the project folder):
#   python3 scripts/make_itol_taxonomy.py [tree] [taxonomy] [output_dir]
#   tree        Newick tree whose tips are portal IDs
#               (default: 06_species_tree/astral/species_tree_astral.nwk).
#   taxonomy    TSV with portal, species, phylum, subphylum, class and order columns
#               (default: 00_inputs/cohort_taxonomy.tsv).
#   output_dir  Where the files are saved (default: 06_species_tree/itol).
#
# Output: itol_<rank>.txt for each rank, and itol_labels.txt, which shows species names instead
# of portal IDs as tip labels. Phylum and subphylum files have a legend; class and order files
# show the group names inside the strip instead, as they have too many groups for a legend.

import csv
import os
import re
import sys

tree_path = sys.argv[1] if len(sys.argv) > 1 else "06_species_tree/astral/species_tree_astral.nwk"
taxonomy_path = sys.argv[2] if len(sys.argv) > 2 else "00_inputs/cohort_taxonomy.tsv"
out_dir = sys.argv[3] if len(sys.argv) > 3 else "06_species_tree/itol"
RANKS = ["phylum", "subphylum", "class", "order"]

# Fixed colours for phyla and subphyla; class and order cycle through a qualitative palette
PHYLUM_COLOURS = {
    "Ascomycota": "#d95f02", "Basidiomycota": "#1b9e77", "Mucoromycota": "#7570b3",
    "Zoopagomycota": "#e7298a", "Blastocladiomycota": "#66a61e", "Chytridiomycota": "#e6ab02",
    "Monoblepharidomycota": "#a6761d", "Neocallimastigomycota": "#666666",
}
SUBPHYLUM_COLOURS = {
    "Pezizomycotina": "#fc8d62", "Saccharomycotina": "#fdbf6f", "Taphrinomycotina": "#b15928",
    "Agaricomycotina": "#66c2a5", "Pucciniomycotina": "#1f78b4", "Ustilaginomycotina": "#a6cee3",
    "Mortierellomycotina": "#cab2d6", "Mucoromycotina": "#6a3d9a", "Glomeromycotina": "#9e9ac8",
    "Entomophthoromycotina": "#e7298a", "none": "#dddddd",
}
CYCLE = ["#8dd3c7", "#fdb462", "#bebada", "#fb8072", "#80b1d3", "#b3de69", "#fccde5",
         "#ffffb3", "#bc80bd", "#ccebc5", "#ffed6f", "#d9d9d9"]


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


for path in (tree_path, taxonomy_path):
    if not os.path.isfile(path):
        fail(f"'{path}' not found\nUsage: {sys.argv[0]} [tree] [taxonomy] [output_dir]")

with open(tree_path) as f:
    newick = re.sub(r"'[^']*'|\[[^\]]*\]", "", f.read())
tips = re.findall(r"[(,]\s*([^,():;\s]+)", newick)

with open(taxonomy_path) as f:
    taxonomy = {row["portal"]: row for row in csv.DictReader(f, delimiter="\t")}
missing = [t for t in tips if t not in taxonomy]
if missing:
    fail(f"{len(missing)} tips are not in '{taxonomy_path}': {', '.join(missing)}")


def colours_for(rank):
    """Colour per group: fixed for phylum and subphylum, cycled in tip order otherwise."""
    fixed = {"phylum": PHYLUM_COLOURS, "subphylum": SUBPHYLUM_COLOURS}.get(rank)
    colours, previous = {}, None
    for tip in tips:
        group = taxonomy[tip][rank]
        if group in colours:
            previous = colours[group]
            continue
        if fixed is not None:
            if group not in fixed:
                fail(f"no colour for {rank} '{group}'; add it to the script")
            colours[group] = fixed[group]
        else:
            colour = CYCLE[len(colours) % len(CYCLE)]
            if colour == previous:
                colour = CYCLE[(len(colours) + 1) % len(CYCLE)]
            colours[group] = colour
        previous = colours[group]
    return colours


os.makedirs(out_dir, exist_ok=True)
for rank in RANKS:
    colours = colours_for(rank)
    lines = ["DATASET_COLORSTRIP", "SEPARATOR TAB", f"DATASET_LABEL\t{rank.capitalize()}",
             "COLOR\t#000000", "STRIP_WIDTH\t30", "MARGIN\t5", "BORDER_WIDTH\t0"]
    if rank in ("phylum", "subphylum"):
        groups = list(colours)
        lines += [f"LEGEND_TITLE\t{rank.capitalize()}",
                  "LEGEND_SHAPES\t" + "\t".join("1" for _ in groups),
                  "LEGEND_COLORS\t" + "\t".join(colours[g] for g in groups),
                  "LEGEND_LABELS\t" + "\t".join(groups)]
    else:
        lines += ["SHOW_STRIP_LABELS\t1", "STRIP_LABEL_POSITION\tcenter", "STRIP_LABEL_SIZE_FACTOR\t0.8",
                  "STRIP_LABEL_COLOR\t#000000", "STRIP_LABEL_ROTATION\t0"]
    lines.append("DATA")
    lines += [f"{tip}\t{colours[taxonomy[tip][rank]]}\t{taxonomy[tip][rank]}" for tip in tips]
    with open(os.path.join(out_dir, f"itol_{rank}.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"  {rank}: {len(colours)} groups -> {out_dir}/itol_{rank}.txt")

with open(os.path.join(out_dir, "itol_labels.txt"), "w") as f:
    f.write("LABELS\nSEPARATOR TAB\nDATA\n")
    f.write("".join(f"{tip}\t{taxonomy[tip]['species']} ({tip})\n" for tip in tips))
print(f"Annotation for {len(tips)} tips saved to {out_dir}/")
