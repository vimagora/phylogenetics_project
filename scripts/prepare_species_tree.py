#!/usr/bin/env python3
# prepare_species_tree.py
# Description: Turns a rooted species tree (by default the ASTRAL tree from species_trees.sbatch)
#              into the plain rooted topology OrthoFinder takes with -s: branch lengths, support
#              values and internal node labels are removed. Checks that the tree is rooted
#              (two children at the root) and that its tips are exactly the cohort portals,
#              which are also OrthoFinder's species names (proteome file names).
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-20
# Version: 1.0
#
# Usage (from the project folder):
#   python3 scripts/prepare_species_tree.py [tree] [portal_list] [output_tree]
#   tree         Rooted Newick tree (default: 06_species_tree/astral/species_tree_astral.nwk).
#   portal_list  TSV with a header row and portal IDs in column 1 (default: 00_inputs/cohort_datasets.tsv).
#   output_tree  Where the topology is saved (default: 07_orthofinder/species_tree_rooted.nwk).

import os
import re
import sys

tree_path = sys.argv[1] if len(sys.argv) > 1 else "06_species_tree/astral/species_tree_astral.nwk"
portal_list = sys.argv[2] if len(sys.argv) > 2 else "00_inputs/cohort_datasets.tsv"
out_path = sys.argv[3] if len(sys.argv) > 3 else "07_orthofinder/species_tree_rooted.nwk"


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


for path in (tree_path, portal_list):
    if not os.path.isfile(path):
        fail(f"'{path}' not found\nUsage: {sys.argv[0]} [tree] [portal_list] [output_tree]")

with open(portal_list) as f:
    portals = {line.split("\t")[0].strip() for line in f.read().splitlines()[1:]} - {""}

with open(tree_path) as f:
    newick = "".join(f.read().split())
if newick.count(";") != 1 or not newick.endswith(";"):
    fail(f"'{tree_path}' should hold exactly one Newick tree")

# Remove quoted or bracketed annotations, branch lengths, and labels after ')' (support values)
topology = re.sub(r"'[^']*'|\[[^\]]*\]", "", newick)
topology = re.sub(r":[^,();]*", "", topology)
topology = re.sub(r"\)[^,();]+", ")", topology)

tips = re.findall(r"[(,]([^,();]+)", topology)
if len(tips) != len(set(tips)):
    fail("some tips appear twice in the tree")
if set(tips) != portals:
    missing, extra = sorted(portals - set(tips)), sorted(set(tips) - portals)
    fail(f"tips do not match the portal list; missing from the tree: {missing}; not on the list: {extra}")

# Rooted means the outermost parentheses hold exactly two children
depth, root_children = 0, 1
for char in topology[1:-2]:
    if char == "(":
        depth += 1
    elif char == ")":
        depth -= 1
    elif char == "," and depth == 0:
        root_children += 1
if not (topology.startswith("(") and topology.endswith(");")) or root_children != 2:
    fail(f"the tree is not rooted ({root_children} children at the root); root it on the outgroup first")

os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
with open(out_path, "w") as f:
    f.write(topology + "\n")
print(f"Rooted topology of {len(tips)} species saved to {out_path}")
