#!/usr/bin/env python3
# structure_trees.py
# Description: Foldtree-style structure trees: an all-vs-all Foldseek comparison of the trimmed
#              models, then one neighbour-joining tree per structural similarity score:
#              fident (sequence identity of the structural alignment), alntmscore (TM-score of the
#              alignment) and lddt. Each pair's score is the mean of both directions, the distance
#              is 1 - score, and trees are rooted at their midpoint.
#              Method after Moi et al. (2025, Nat Struct Mol Biol, "Foldtree").
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-22
# Version: 1.0
#
# Usage (from the project folder; needs foldseek and Biopython):
#   python3 scripts/structure_trees.py [--models-dir DIR] [--out-dir DIR] [--threads N]
#   --models-dir  Trimmed models from fetch_afdb_models.py (default: 08_structures/models_trimmed).
#   --out-dir     Where results are saved (default: 08_structures/trees).
#
# Output, inside out-dir: foldseek_allvsall.tsv, and structure_<score>.nwk for each score.

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

SCORES = ["fident", "alntmscore", "lddt"]

parser = argparse.ArgumentParser(description="Foldtree-style structure trees")
parser.add_argument("--models-dir", default="08_structures/models_trimmed")
parser.add_argument("--out-dir", default="08_structures/trees")
parser.add_argument("--threads", type=int, default=os.cpu_count())
args = parser.parse_args()

if shutil.which("foldseek") is None:
    sys.exit("ERROR: foldseek not found; activate the structure_test environment")
models = sorted(f[:-4] for f in os.listdir(args.models_dir) if f.endswith(".pdb"))
if len(models) < 4:
    sys.exit(f"ERROR: fewer than 4 models in '{args.models_dir}', run scripts/fetch_afdb_models.py first")
os.makedirs(args.out_dir, exist_ok=True)

# All-vs-all search: exhaustive and without an E-value cut-off, so that every pair gets scores
table = os.path.join(args.out_dir, "foldseek_allvsall.tsv")
with tempfile.TemporaryDirectory() as tmp:
    subprocess.run(["foldseek", "easy-search", args.models_dir, args.models_dir, table, tmp,
                    "--exhaustive-search", "1", "-e", "inf", "--threads", str(args.threads),
                    "--format-output", "query,target," + ",".join(SCORES)],
                   check=True, stdout=subprocess.DEVNULL)


def model_name(name):
    """Foldseek names entries after the file, sometimes with the extension or a chain suffix."""
    for suffix in (".pdb_A", "_A", ".pdb"):
        if name.endswith(suffix) and name[:-len(suffix)] in index:
            return name[:-len(suffix)]
    return name


index = {m: i for i, m in enumerate(models)}
scores = {s: {} for s in SCORES}
with open(table) as f:
    for line in f:
        query, target, *values = line.rstrip("\n").split("\t")
        query, target = model_name(query), model_name(target)
        if query not in index or target not in index:
            sys.exit(f"ERROR: unexpected entry name in Foldseek output: {query} / {target}")
        for score, value in zip(SCORES, values):
            pair = scores[score].setdefault(frozenset((query, target)), [])
            pair.append(float(value))

n_pairs = len(models) * (len(models) - 1) // 2
for score in SCORES:
    found = sum(1 for pair in scores[score] if len(pair) == 2)
    if found < n_pairs:
        print(f"  note: {n_pairs - found} of {n_pairs} pairs have no {score} hit; their distance is set to 1")
    # Lower-triangular distance matrix, 1 - mean score of both directions (0 on the diagonal)
    matrix = []
    for i, a in enumerate(models):
        row = []
        for b in models[:i]:
            values = scores[score].get(frozenset((a, b)), [])
            row.append(1 - (sum(values) / len(values)) if values else 1.0)
        row.append(0.0)
        matrix.append(row)
    tree = DistanceTreeConstructor().nj(DistanceMatrix(models, matrix))
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = None
        if clade.branch_length is not None and clade.branch_length < 0:
            clade.branch_length = 0.0
    tree.root_at_midpoint()
    path = os.path.join(args.out_dir, f"structure_{score}.nwk")
    Phylo.write(tree, path, "newick")
    print(f"  {score}: neighbour-joining tree of {len(models)} models -> {path}")
