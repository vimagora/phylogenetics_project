#!/usr/bin/env python3
# check_subfamily.py
# Description: Checks a subfamily's tree (step F14), for a person to review; nothing is decided here.
#              It reports: whether the ingroup is monophyletic; whether each outgroup lineage is
#              monophyletic; whether rooting on each lineage gives the same ingroup root; whether the
#              outgroup root is accepted: when the ingroup is a clade, so that every outgroup tip falls
#              outside it, even if a lineage is not a clade itself (that is reported); where the midpoint and
#              MAD roots fall against the outgroup root, and whether they only move the outgroups to the other
#              side of the ingroup's first split; ingroup support; and how many ingroup
#              splits the whole-family trees share. MAD (minimal ancestor deviation; Tria, Landan & Dagan
#              2017, Nat Ecol Evol 1:0193) assumes roughly clock-like evolution, so a disagreement is
#              something to explain, not a reason to replace the outgroup root.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-11-05
# Version: 1.1
#
# Usage (from the project folder; needs numpy, e.g. module load python-data):
#   python3 scripts/check_subfamily.py <family> <subfamily> [--tree fasttree|ml]
#   --tree  The subfamily tree to check (default: fasttree).
#   MAD compares all pairs of tips on every branch: seconds for hundreds of tips, minutes for thousands.
#
# Reads, inside <OUT_DIR>/subfamilies/<subfamily> (OUT_DIR from family.conf): tips.tsv, trimmed.faa, and
# fasttree.nwk or ml.treefile; and <OUT_DIR>/tree/fasttree.nwk and ml.treefile, if present.
#
# Output, in the same folder: checks_<tree>.md, also printed; and <tree>_rooted.nwk, rooted on the
# outgroups, with the original supports, if the ingroup is monophyletic.

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_config import fail, read_config, read_fasta, read_tsv  # noqa: E402
from newick_tree import read_tree  # noqa: E402

TREES = {"fasttree": "fasttree.nwk", "ml": "ml.treefile"}
HIGH_SUPPORT = {"fasttree": ("≥ 0.9", lambda label: float(label) >= 0.9),
                "ml": ("SH-aLRT ≥ 80 and UFBoot ≥ 95",
                       lambda label: float(label.split("/")[0]) >= 80 and float(label.split("/")[1]) >= 95)}
SAME_POINT = 1e-12


def mad_roots(tree):
    """All root positions as (deviation, node, distance above the node, point), best first. A point names
    where the root falls, since the branches around a node tie when the best position is the node itself."""
    nodes = tree.preorder()
    index = {id(n): i for i, n in enumerate(nodes)}
    tips = list(tree.tips)
    # distance from every tip to every node: first to the nodes above the tip, then to all the others
    under = np.zeros((len(tips), len(nodes)), bool)
    dist = np.zeros((len(tips), len(nodes)))
    for i, tip in enumerate(tips):
        under[i, index[id(tree.tips[tip])]] = True
    for node in reversed(nodes):
        if node.parent:
            n, p = index[id(node)], index[id(node.parent)]
            dist[under[:, n], p] = dist[under[:, n], n] + node.length
            under[:, p] |= under[:, n]
    for node in nodes[1:]:
        n = index[id(node)]
        dist[~under[:, n], n] = dist[~under[:, n], index[id(node.parent)]] + node.length
    b, c = np.triu_indices(len(tips), 1)
    tip_dist = dist[:, [index[id(tree.tips[t])] for t in tips]][b, c]
    tip_dist = np.where(tip_dist > 0, tip_dist, 1e-12)
    roots = []
    for node in nodes[1:]:
        if node.length <= 0:
            continue
        n, p = index[id(node)], index[id(node.parent)]
        below_b, below_c = under[b, n], under[c, n]
        across = below_b != below_c
        d = np.where(below_b & below_c, dist[b, n] - dist[c, n], dist[b, p] - dist[c, p]) / tip_dist
        to_node = np.where(below_b, dist[b, n], dist[c, n])[across]
        across_dist = tip_dist[across]
        rho = -np.sum((2 * to_node - across_dist) / across_dist ** 2) / (2 * np.sum(1 / across_dist ** 2))
        rho = float(np.clip(rho, 0, node.length))
        d[across] = (2 * (to_node + rho) - across_dist) / across_dist
        if rho <= SAME_POINT:
            point = ("node", n)
        elif rho >= node.length - SAME_POINT:
            point = ("node", p)
        else:
            point = ("branch", n)
        roots.append((float(np.sqrt(np.mean(d ** 2))), node, rho, point))
    roots.sort(key=lambda r: r[0])
    return roots


def splits(tree, keep):
    """The non-trivial splits of a tree pruned to the tips in keep, each as the integer whose bits are the
    tips of the side without the first tip."""
    bit = {t: 1 << i for i, t in enumerate(sorted(keep))}
    full, found, masks = (1 << len(keep)) - 1, set(), {}
    for node in reversed(tree.preorder()):
        mask = bit.get(node.name, 0)
        for child in node.children:
            mask |= masks.pop(id(child))
        masks[id(node)] = mask
        if mask & 1:
            mask ^= full
        if bin(mask).count("1") >= 2 and bin(full ^ mask).count("1") >= 2:
            found.add(mask)
    return found


parser = argparse.ArgumentParser(description="Checks of a subfamily's tree")
parser.add_argument("family")
parser.add_argument("subfamily")
parser.add_argument("--tree", choices=TREES, default="fasttree")
args = parser.parse_args()
config = read_config(args.family)
out_dir = config["OUT_DIR"]
sub_dir = os.path.join(out_dir, "subfamilies", args.subfamily)
tree_path = os.path.join(sub_dir, TREES[args.tree])
if not os.path.isfile(tree_path):
    fail(f"'{tree_path}' not found")

tips_table = {t["tip_id"]: t for t in read_tsv(os.path.join(sub_dir, "tips.tsv"))}
tree = read_tree(tree_path)
if set(tree.tips) != set(tips_table):
    fail(f"the tips of {tree_path} differ from {sub_dir}/tips.tsv; rebuild the tree")
ingroup = [t for t, info in tips_table.items() if info["set"] == "ingroup"]
ingroup_set = set(ingroup)
outgroups = [t for t in tips_table if t not in ingroup_set]
if not outgroups:
    fail(f"{args.subfamily} has no outgroups in tips.tsv")
lineages = {}
for tip in outgroups:
    lineages.setdefault(tips_table[tip]["set"], []).append(tip)
references = {r["reference_id"]: r for r in read_tsv(os.path.join(config["FAMILY_DIR"], "references.tsv"))}
report = [f"# Checks of the {args.subfamily} tree ({args.tree})", "",
          f"{len(ingroup)} ingroup tips; outgroups: " + ", ".join(f"{k} {len(v)}" for k, v in lineages.items()) + "."]
trimmed = os.path.join(sub_dir, "trimmed.faa")
if os.path.isfile(trimmed):
    report.append(f"Alignment: {len(next(iter(read_fasta(trimmed).values())))} columns after trimming.")


def named_in(tip_list):
    """The orthogroups of some tips and the names of the references among them, for reading a root."""
    groups = {}
    for tip in tip_list:
        group = tips_table[tip]["orthogroup"] or "(none)"
        groups[group] = groups.get(group, 0) + 1
    main = sorted(groups, key=lambda g: (-groups[g], g))[:3]
    names = sorted({references[r]["protein_name"] for t in tip_list for r in tips_table[t]["references"].split(",")
                    if r in references and references[r]["evidence"] == "characterized"})
    return ", ".join(f"{g} {groups[g]}" for g in main) + (f"; names: {', '.join(names[:12])}" if names else "")


# 1. The ingroup, and the same ingroup root from each outgroup lineage
report += ["", "## Monophyly and root", ""]
ingroup_roots = {}
ingroup_clade = True
for lineage, members in lineages.items():
    rooted = read_tree(tree_path)
    rooted.root_on_tip(members[0])
    node = rooted.mrca(ingroup)
    inside = rooted.tips_under(node)
    intruders = {}
    for tip in inside:
        if tip not in ingroup_set:
            intruders[tips_table[tip]["set"]] = intruders.get(tips_table[tip]["set"], 0) + 1
    if intruders:
        ingroup_clade = False
    ingroup_roots[lineage] = frozenset(frozenset(ingroup_set.intersection(rooted.tips_under(c)))
                                       for c in node.children) - {frozenset()}
    state = "PASS" if not intruders else "FAIL: inside it, " + ", ".join(f"{k} {v}" for k, v in intruders.items())
    report.append(f"- Ingroup monophyletic, seen from {lineage}: {state} (support {node.label or 'none'})")
same = len(set(ingroup_roots.values())) == 1
report.append(f"- Same ingroup root from all lineages: {'PASS' if same else 'FAIL'}"
              + (" (only one lineage)" if len(lineages) == 1 else ""))

# 2. Each outgroup lineage, then the verdict on the outgroup root
failed_lineages = []
for lineage, members in lineages.items():
    rooted = read_tree(tree_path)
    rooted.root_on_tip(members[0])
    rest = [t for t in tips_table if t not in set(members)]
    node = rooted.mrca(rest)
    clade = len(rooted.tips_under(node)) == len(rest)
    if not clade:
        failed_lineages.append(lineage)
    report.append(f"- Lineage {lineage} monophyletic: {'PASS' if clade else 'FAIL'}"
                  + (f" (support {node.label or 'none'})" if clade else ""))

if ingroup_clade:
    report.append(f"- Outgroup root: ACCEPTED, all {len(outgroups)} outgroup tips fall outside the ingroup"
                  + (f" (not a clade itself, reported only: {', '.join(failed_lineages)})" if failed_lineages else ""))
else:
    report.append("- Outgroup root: NOT ACCEPTED, outgroup tips fall inside the ingroup")

# 3. The midpoint and MAD roots against the outgroup root, and the ingroup's first split under that root
report += ["", "## Other roots", ""]
first_split = set()
if ingroup_clade:
    rooted = read_tree(tree_path)
    rooted.root_on_clade(outgroups)
    first_split = {frozenset(rooted.tips_under(c)) & ingroup_set for c in rooted.mrca(ingroup).children}


def place_root(method, side, extra=""):
    side = set(side)
    one = [t for t in ingroup if t in side]
    other = [t for t in ingroup if t not in side]
    if not one or not other:
        report.append(f"- {method}: outside the ingroup, as the outgroup root{extra}")
        return
    smaller, larger = sorted((one, other), key=len)
    same_split = frozenset(smaller) in first_split or frozenset(larger) in first_split
    report.append(f"- {method}: INSIDE the ingroup, separating {len(smaller)} tip(s) from {len(larger)}{extra}. "
                  f"The {len(smaller)}: {named_in(smaller)}"
                  + (". This is the ingroup's first split under the outgroup root: the two roots agree on it and "
                     "differ only in the side the outgroups join" if same_split else ""))


midpoint = read_tree(tree_path)
midpoint.midpoint_root()
place_root("Midpoint root", midpoint.tips_under(midpoint.root.children[0]))
unrooted = read_tree(tree_path)
roots = mad_roots(unrooted)
best = roots[0]
second = next((r[0] for r in roots[1:] if r[3] != best[3]), None)
ambiguity = f"; ambiguity index {best[0] / second:.3f}, where 1 means unresolved" if second else ""
place_root("MAD root", unrooted.tips_under(best[1]), f" (deviation {best[0]:.4f}{ambiguity})")

# 4. Support within the ingroup, and the rooted tree
rooted = read_tree(tree_path)
if ingroup_clade and rooted.root_on_clade(outgroups):
    rooted_path = os.path.join(sub_dir, f"{args.tree}_rooted.nwk")
    with open(rooted_path, "w") as f:
        f.write(rooted.newick() + "\n")
    report += ["", "## Support", ""]
    labelled, high = 0, 0
    stack = [rooted.mrca(ingroup)]
    text, is_high = HIGH_SUPPORT[args.tree]
    while stack:
        node = stack.pop()
        stack.extend(c for c in node.children if c.children)
        if node.label:
            labelled += 1
            high += is_high(node.label)
    report.append(f"- {high} of {labelled} ingroup nodes ({high / max(labelled, 1):.0%}) have support {text}")
    report.append(f"- Rooted tree: {rooted_path}")
else:
    report += ["", "No rooted tree was written, since the ingroup is not monophyletic."]

# 5. The same tips in the whole-family trees, as a check, not a requirement
report += ["", "## Against the whole-family trees (ingroup tips only)", ""]
own = splits(read_tree(tree_path), ingroup_set)
compared = False
for name, file_name in TREES.items():
    path = os.path.join(out_dir, "tree", file_name)
    if os.path.isfile(path):
        family_tree = read_tree(path)
        missing = ingroup_set - set(family_tree.tips)
        if missing:
            report.append(f"- {name}: {len(missing)} ingroup tips are not in this tree; not compared")
        else:
            shared = len(own & splits(family_tree, ingroup_set))
            report.append(f"- {name}: {shared} of {len(own)} splits shared ({shared / max(len(own), 1):.0%})")
        compared = True
if not compared:
    report.append("- No whole-family tree found")

path = os.path.join(sub_dir, f"checks_{args.tree}.md")
with open(path, "w", encoding="utf-8") as f:
    f.write("\n".join(report) + "\n")
print("\n".join(report))
print(f"\nWritten to {path}")
