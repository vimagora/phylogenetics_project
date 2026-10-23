#!/usr/bin/env python3
# prepare_species_tree.py
# Description: Turns a rooted species tree (by default the ASTRAL tree from species_trees.sbatch)
#              into the plain rooted topology OrthoFinder takes with -s: branch lengths, support
#              values and internal node labels are removed. Checks that the tree is rooted
#              (two children at the root) and that its tips are exactly the cohort portals,
#              which are also OrthoFinder's species names (proteome file names). With --outgroup,
#              the tree is first rerooted on the branch that separates the outgroup from the rest.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-20
# Version: 1.1
#
# Usage (from the project folder):
#   python3 scripts/prepare_species_tree.py [tree] [portal_list] [output_tree] [--outgroup PORTAL,...]
#   tree         Newick tree, rooted unless --outgroup is given
#                (default: 06_species_tree/astral/species_tree_astral.nwk).
#   portal_list  TSV with a header row and portal IDs in column 1 (default: 00_inputs/cohort_datasets.tsv).
#   output_tree  Where the topology is saved (default: 07_orthofinder/species_tree_rooted.nwk).
#   --outgroup   Portal(s) to root on, comma-separated; they must form a clade in the unrooted tree.
#                Without it, the tree's own root is kept. For OrthoFinder we root on Catan2
#                (Blastocladiomycota diverging first, as in the concatenated tree of Li et al. 2021);
#                the species trees themselves are rooted on Spipu1, which falls inside Chytridiomyceta.

import argparse
import os
import re
import sys

parser = argparse.ArgumentParser(description="Rooted species-tree topology for OrthoFinder")
parser.add_argument("tree", nargs="?", default="06_species_tree/astral/species_tree_astral.nwk")
parser.add_argument("portal_list", nargs="?", default="00_inputs/cohort_datasets.tsv")
parser.add_argument("output_tree", nargs="?", default="07_orthofinder/species_tree_rooted.nwk")
parser.add_argument("--outgroup", default="")
args = parser.parse_args()
tree_path, portal_list, out_path = args.tree, args.portal_list, args.output_tree


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


for path in (tree_path, portal_list):
    if not os.path.isfile(path):
        fail(f"'{path}' not found\nUsage: {sys.argv[0]} [tree] [portal_list] [output_tree] [--outgroup PORTAL,...]")

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



def parse(text):
    """Nested lists of tip names from a plain Newick topology."""
    stack, current, name = [], [], ""
    for char in text.rstrip(";"):
        if char == "(":
            stack.append(current)
            current = []
        elif char in ",)":
            if name:
                current.append(name)
                name = ""
            if char == ")":
                node, current = current, stack.pop()
                current.append(node)
        else:
            name += char
    return current[0]


def to_newick(node):
    return node if isinstance(node, str) else "(" + ",".join(to_newick(child) for child in node) + ")"


def tips_of(node):
    return {node} if isinstance(node, str) else set().union(*(tips_of(child) for child in node))


def reroot(tree, outgroup):
    """Tree rooted on the branch above the outgroup clade, or None if the outgroup is not a clade."""
    # Unroot first: a root with two children becomes one node holding the contents of one of them
    if len(tree) == 2:
        if outgroup in (tips_of(tree[0]), tips_of(tree[1])):
            return tree
        inner, other = (tree[0], tree[1]) if not isinstance(tree[0], str) else (tree[1], tree[0])
        tree = [*inner, other]

    # Walk down the tree; when the root moves above a child, the rest of the tree becomes its sister
    def search(node, above):
        for child in node:
            rest = [c for c in node if c is not child] + ([above] if above else [])
            if tips_of(child) == outgroup:
                return [child, rest]
            if not isinstance(child, str):
                found = search(child, rest)
                if found:
                    return found
        return None

    return search(tree, None)


if args.outgroup:
    outgroup = {portal.strip() for portal in args.outgroup.split(",") if portal.strip()}
    if not outgroup <= set(tips):
        fail(f"outgroup portals not in the tree: {sorted(outgroup - set(tips))}")
    tree = parse(topology)
    rooted = reroot(tree, outgroup) or reroot(tree, set(tips) - outgroup)
    if rooted is None:
        fail(f"the outgroup {sorted(outgroup)} is not a clade in the tree")
    if tips_of(rooted[0]) != outgroup:
        rooted = rooted[::-1]
    topology = to_newick(rooted) + ";"
    if sorted(re.findall(r"[(,]([^,();]+)", topology)) != sorted(tips):
        fail("rerooting changed the tips; this is a bug")

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
rooted_on = f", rooted on {args.outgroup}" if args.outgroup else ""
print(f"Rooted topology of {len(tips)} species{rooted_on} saved to {out_path}")
