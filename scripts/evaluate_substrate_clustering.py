#!/usr/bin/env python3
# evaluate_substrate_clustering.py
# Description: Scores how well trees of sugar transporters group proteins with the same substrates,
#              using only characterized transporters (sensors, D. squalens and transporters without
#              a specific substrate are left out of the scores):
#              (a) nearest neighbour: share of transporters whose nearest characterized neighbour
#                  (shortest path through the tree) has at least one substrate in common, with a
#                  p-value from shuffling the substrate labels over the tips;
#              (b) purity: for each substrate named for >= 3 transporters, the smallest clade (tree
#                  rooted at its midpoint) holding all of them, and the share of the clade's
#                  characterized transporters that take up that substrate; with the clade's support
#                  where the tree has support values.
#              Generic labels stand for all their members when matching: hexose = D-glucose,
#              D-fructose, D-mannose, D-galactose; pentose = D-xylose, L-arabinose; polyol = glycerol
#              and the sugar alcohols; cellodextrin = cellobiose, cellotriose, cellotetraose.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-22
# Version: 1.0
#
# Usage (from the project folder; needs Biopython):
#   python3 scripts/evaluate_substrate_clustering.py [tree.nwk ...] [--out-dir DIR] [--permutations 1000]
#   Trees default to the four in 08_structures/trees (structure_fident, structure_alntmscore,
#   structure_lddt, sequence_ml.treefile).
#
# Output, inside out-dir (default: 08_structures/evaluation):
#   summary.tsv              per tree: (a) share, mean share with shuffled labels and p-value, (b) mean purity
#   per_substrate.tsv        per tree and substrate: transporters, clade size, purity, support
#   nearest_neighbours.tsv   per tree and transporter: nearest neighbour and whether they share a substrate
#   itol_substrates.txt      iTOL binary dataset: which substrates each protein takes up
#   itol_category.txt        iTOL colour strip: characterized transporter, sensor, D. squalens, unspecified
#   itol_labels.txt          iTOL tip labels: protein name and species

import argparse
import csv
import os
import random
import re
import sys

from Bio import Phylo

IDS_TABLE = "08_structures/structure_ids.tsv"
DEFAULT_TREES = ["08_structures/trees/structure_fident.nwk", "08_structures/trees/structure_alntmscore.nwk",
                 "08_structures/trees/structure_lddt.nwk", "08_structures/trees/sequence_ml.treefile"]
# Spelling differences, including the shorter names of the predicted D. squalens specificities
SYNONYMS = {"galactose": "d-galactose", "mannitol": "d-mannitol", "inositol": "myo-inositol",
            "xylose": "d-xylose", "arabitol": "d-arabitol", "uronic": "d-galacturonic acid"}
GENERIC = {
    "hexose": {"d-glucose", "d-fructose", "d-mannose", "d-galactose"},
    "pentose": {"d-xylose", "l-arabinose"},
    "polyol": {"glycerol", "d-sorbitol", "d-mannitol", "d-arabitol", "l-arabitol", "galactitol", "xylitol",
               "ribitol"},
    "cellodextrin": {"cellobiose", "cellotriose", "cellotetraose"},
}
NOT_SPECIFIC = {"sugar (unspecified)", ""}
# Columns of the iTOL substrate dataset: (label, substrates that set it), with a colour each
ITOL_COLUMNS = [
    ("D-glucose", {"d-glucose"}, "#1f77b4"), ("D-fructose", {"d-fructose"}, "#aec7e8"),
    ("D-mannose", {"d-mannose"}, "#17becf"), ("D-galactose", {"d-galactose"}, "#9edae5"),
    ("D-xylose", {"d-xylose"}, "#2ca02c"), ("L-arabinose", {"l-arabinose"}, "#98df8a"),
    ("L-rhamnose", {"l-rhamnose"}, "#bcbd22"),
    ("cello-/xylo-oligosaccharides", {"cellobiose", "cellotriose", "cellotetraose", "sophorose", "cellobionic acid",
                                     "xylobiose"}, "#ff7f0e"),
    ("lactose/galactosides", {"lactose", "galactoside"}, "#ffbb78"),
    ("maltose/alpha-glucosides", {"maltose", "alpha-glucoside"}, "#d62728"), ("sucrose", {"sucrose"}, "#ff9896"),
    ("inositols", {"myo-inositol", "d-chiro-inositol"}, "#9467bd"), ("glycerol", {"glycerol"}, "#c5b0d5"),
    ("other polyols", GENERIC["polyol"] - {"glycerol"}, "#8c564b"),
    ("uronic acids", {"d-galacturonic acid", "d-glucuronic acid", "alginate"}, "#e377c2"),
    ("quinic acid", {"quinic acid"}, "#7f7f7f"),
]


def tokens(substrates):
    """Normalized substrate names of a ';'- or ','-separated label."""
    names = {SYNONYMS.get(t.strip().lower(), t.strip().lower()) for t in re.split(r"[;,]", substrates)}
    return names - NOT_SPECIFIC


def expanded(names):
    """Specific substrates: generic labels replaced by their members."""
    out = set()
    for name in names:
        out |= GENERIC.get(name, {name})
    return out


def path_lengths(tree, tips):
    """All pairwise path lengths between the given tips."""
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            parents[child] = clade
    depth = tree.depths(unit_branch_lengths=False)
    paths = {}
    for tip in tips:
        node, path = tip, []
        while node is not None:
            path.append(node)
            node = parents.get(node)
        paths[tip] = path
    ancestors = {tip: set(map(id, path)) for tip, path in paths.items()}
    distances = {}
    for i, a in enumerate(tips):
        for b in tips[i + 1:]:
            lca = next(node for node in paths[b] if id(node) in ancestors[a])
            distances[(a, b)] = distances[(b, a)] = depth[a] + depth[b] - 2 * depth[lca]
    return distances


def nearest_neighbour_share(tips, labels, distances):
    """Share of tips whose nearest tip(s) share at least one substrate, and the neighbours."""
    hits, details = 0, []
    for a in tips:
        others = [b for b in tips if b is not a]
        best = min(distances[(a, b)] for b in others)
        nearest = [b for b in others if distances[(a, b)] <= best + 1e-9]
        shared = any(labels[a] & labels[b] for b in nearest)
        hits += shared
        details.append((a, nearest, shared))
    return hits / len(tips), details


parser = argparse.ArgumentParser(description="Score substrate clustering on transporter trees")
parser.add_argument("trees", nargs="*", default=DEFAULT_TREES)
parser.add_argument("--out-dir", default="08_structures/evaluation")
parser.add_argument("--permutations", type=int, default=1000)
args = parser.parse_args()

with open(IDS_TABLE) as f:
    proteins = {row["structure_id"]: row for row in csv.DictReader(f, delimiter="\t")}


def category(row):
    if row["source"] == "dicsqu":
        return "D. squalens (predicted)"
    if row["role"] == "sensor":
        return "sensor"
    if not tokens(row["substrates"]):
        return "substrate unspecified"
    return "characterized transporter"


os.makedirs(args.out_dir, exist_ok=True)
summary, per_substrate, neighbours = [], [], []
random.seed(1)
for tree_path in args.trees:
    if not os.path.isfile(tree_path):
        sys.exit(f"ERROR: tree '{tree_path}' not found")
    tree_name = os.path.basename(tree_path).rsplit(".", 1)[0]
    tree = Phylo.read(tree_path, "newick")
    tree.root_at_midpoint()
    unknown = [t.name for t in tree.get_terminals() if t.name not in proteins]
    if unknown:
        sys.exit(f"ERROR: {tree_name} has tips not in {IDS_TABLE}: {', '.join(unknown[:5])}")
    tips = [t for t in tree.get_terminals() if category(proteins[t.name]) == "characterized transporter"]
    named = {t: tokens(proteins[t.name]["substrates"]) for t in tips}
    labels = {t: expanded(named[t]) for t in tips}
    distances = path_lengths(tree, tips)

    # (a) nearest neighbour, against shuffled labels
    share, details = nearest_neighbour_share(tips, labels, distances)
    values = list(labels.values())
    shuffled = []
    for _ in range(args.permutations):
        random.shuffle(values)
        shuffled.append(nearest_neighbour_share(tips, dict(zip(tips, values)), distances)[0])
    p_value = (sum(x >= share for x in shuffled) + 1) / (args.permutations + 1)
    chance = sum(shuffled) / len(shuffled)
    for a, nearest, shared in details:
        neighbours.append((tree_name, a.name, ",".join(b.name for b in nearest), shared))

    # (b) purity of the smallest clade holding all transporters named for a substrate
    purities = []
    for substrate in sorted({s for names in named.values() for s in names} - set(GENERIC)):
        carriers = [t for t in tips if substrate in named[t]]
        if len(carriers) < 3:
            continue
        clade = tree.common_ancestor(carriers)
        members = [t for t in clade.get_terminals() if t in labels]
        purity = sum(substrate in labels[t] for t in members) / len(members)
        support = clade.confidence if clade.confidence is not None and clade is not tree.root else ""
        purities.append(purity)
        per_substrate.append((tree_name, substrate, len(carriers), len(members), round(purity, 3), support))

    mean_purity = sum(purities) / len(purities)
    summary.append((tree_name, len(tips), round(share, 3), round(chance, 3), round(p_value, 4), len(purities),
                    round(mean_purity, 3)))
    print(f"  {tree_name}: (a) nearest neighbour shares a substrate for {share:.1%} of {len(tips)} "
          f"transporters ({chance:.1%} with shuffled labels, p = {p_value:.4f}); "
          f"(b) mean purity {mean_purity:.2f} over {len(purities)} substrates")


def write(name, header, rows):
    with open(os.path.join(args.out_dir, name), "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


write("summary.tsv", ["tree", "transporters_scored", "nn_share", "nn_share_shuffled", "nn_p_value", "substrates_scored", "mean_purity"],
      summary)
write("per_substrate.tsv", ["tree", "substrate", "transporters", "clade_transporters", "purity", "clade_support"],
      per_substrate)
write("nearest_neighbours.tsv", ["tree", "structure_id", "nearest", "shares_substrate"], neighbours)

# iTOL files for all proteins in the trees
in_trees = sorted({t.name for path in args.trees for t in Phylo.read(path, "newick").get_terminals()})
lines = ["DATASET_BINARY", "SEPARATOR TAB", "DATASET_LABEL\tSubstrates", "COLOR\t#000000",
         "FIELD_SHAPES\t" + "\t".join("1" for _ in ITOL_COLUMNS),
         "FIELD_LABELS\t" + "\t".join(c[0] for c in ITOL_COLUMNS),
         "FIELD_COLORS\t" + "\t".join(c[2] for c in ITOL_COLUMNS), "SHOW_LABELS\t1", "DATA"]
for sid in in_trees:
    names = expanded(tokens(proteins[sid]["substrates"].removeprefix("predicted: ")))
    lines.append(sid + "\t" + "\t".join("1" if names & c[1] else "-1" for c in ITOL_COLUMNS))
with open(os.path.join(args.out_dir, "itol_substrates.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
colours = {"characterized transporter": "#4daf4a", "sensor": "#e41a1c", "D. squalens (predicted)": "#377eb8",
           "substrate unspecified": "#999999"}
lines = ["DATASET_COLORSTRIP", "SEPARATOR TAB", "DATASET_LABEL\tCategory", "COLOR\t#000000",
         "LEGEND_TITLE\tCategory", "LEGEND_SHAPES\t" + "\t".join("1" for _ in colours),
         "LEGEND_COLORS\t" + "\t".join(colours.values()), "LEGEND_LABELS\t" + "\t".join(colours), "DATA"]
lines += [f"{sid}\t{colours[category(proteins[sid])]}\t{category(proteins[sid])}" for sid in in_trees]
with open(os.path.join(args.out_dir, "itol_category.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
with open(os.path.join(args.out_dir, "itol_labels.txt"), "w") as f:
    f.write("LABELS\nSEPARATOR TAB\nDATA\n")
    for sid in in_trees:
        species = proteins[sid]["species"].split(" (")[0]
        f.write(f"{sid}\t{proteins[sid]['protein_name']} ({species})\n")
print(f"Results in {args.out_dir}/")
