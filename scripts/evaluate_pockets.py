#!/usr/bin/env python3
# evaluate_pockets.py
# Family-specific: sugar transporters (structure concept test; see README, "Sugar transporter structures").
# Description: Scores whether sugar-binding pocket residues (map_pockets.py) group characterized
#              sugar transporters by substrate better than the whole protein does. Each method gives
#              a distance between every two transporters:
#                pocket_core, pocket_extended  BLOSUM62 distance between pocket signatures
#                all_columns                   the same distance over every well-aligned column
#                sequence_ml                   path length in the maximum-likelihood sequence tree
#              Scores, over the characterized transporters with a model (as in
#              evaluate_substrate_clustering.py):
#                (a) nearest neighbour shares a substrate (ties share the credit), with shuffled labels
#                (a') the same, among transporters < max-identity identical to it (other lineages)
#                (c) the most common substrate of the 3 nearest neighbours is one of its own
#              Baselines: the sequence tree (paired sign-flip test per transporter) and random sets of
#              well-aligned columns of the pocket's size.
#              Decision rule, set before the run: the core pocket is worth pursuing if, on (a) and on
#              (a'), it scores higher than the sequence tree (paired p < 0.05) and higher than 95% of
#              the random column sets.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-23
# Version: 1.0
#
# Usage (from the project folder; needs numpy and Biopython):
#   python3 scripts/evaluate_pockets.py [--max-identity 0.4] [--random-sets 1000] [--permutations 1000]
#
# Output, inside 08_structures/pockets:
#   trees/pocket_core.nwk, trees/pocket_extended.nwk   neighbour-joining trees of all models (midpoint
#                                                       rooted), scored by evaluate_substrate_clustering.py
#                                                       into evaluation/trees/
#   evaluation/summary.tsv          per method: scores, chance, p-values against shuffled labels, the
#                                   sequence tree and the random column sets
#   evaluation/random_columns.tsv   scores of every random column set
#   evaluation/per_transporter.tsv  per method and transporter: credit for (a), (a') and (c), and neighbours

import argparse
import csv
import os
import random
import subprocess
import sys

import numpy as np
from Bio import Phylo
from Bio.Align import substitution_matrices
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from evaluate_substrate_clustering import category, expanded, path_lengths, tokens  # noqa: E402

IDS_TABLE = "08_structures/structure_ids.tsv"
POCKET_DIR = "08_structures/pockets"
SEQUENCE_TREE = "08_structures/trees/sequence_ml.treefile"
EVAL_DIR = os.path.join(POCKET_DIR, "evaluation")
TREE_DIR = os.path.join(POCKET_DIR, "trees")

parser = argparse.ArgumentParser(description="Score pocket signatures against whole-protein baselines")
parser.add_argument("--max-identity", type=float, default=0.4)
parser.add_argument("--random-sets", type=int, default=1000)
parser.add_argument("--permutations", type=int, default=1000)
parser.add_argument("--min-shared", type=float, default=0.5,
                    help="pairs sharing fewer than this share of a signature's columns get no distance")
args = parser.parse_args()

for path in (IDS_TABLE, SEQUENCE_TREE, os.path.join(POCKET_DIR, "pocket_signatures.tsv")):
    if not os.path.isfile(path):
        sys.exit(f"ERROR: '{path}' not found; run scripts/sequence_tree.sh and scripts/map_pockets.py first")

with open(IDS_TABLE) as f:
    proteins = {row["structure_id"]: row for row in csv.DictReader(f, delimiter="\t")}
with open(os.path.join(POCKET_DIR, "pocket_signatures.tsv")) as f:
    signatures = {row["structure_id"]: row for row in csv.DictReader(f, delimiter="\t")}
aligned, name = {}, None
with open(os.path.join(POCKET_DIR, "alignment.faa")) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            aligned[name] = []
        elif name:
            aligned[name].append(line)
aligned = {k: "".join(v) for k, v in aligned.items()}
with open(os.path.join(POCKET_DIR, "core_columns.txt")) as f:
    core_columns = [int(x) - 1 for x in f.read().split()]
with open(os.path.join(POCKET_DIR, "pocket_columns.tsv")) as f:
    pocket_columns = {int(r["column"]) - 1 for r in csv.DictReader(f, delimiter="\t") if r["kept"] == "True"}

models = [sid for sid, row in signatures.items() if row["has_model"] == "True"]
scored = [sid for sid in models if category(proteins[sid]) == "characterized transporter"]
labels = [expanded(tokens(proteins[sid]["substrates"])) for sid in scored]
n = len(scored)

# Residues as indices into BLOSUM62; anything else (gap, X, B, Z) counts as missing
blosum = substitution_matrices.load("BLOSUM62")
letters = "ARNDCQEGHILKMFPSTWYV"
B = np.array([[blosum[a][b] for b in letters] for a in letters], dtype=float)
code = {a: i for i, a in enumerate(letters)}


def encode(strings):
    return np.array([[code.get(c, -1) for c in s] for s in strings], dtype=int)


def blosum_distance(x):
    """1 - normalized BLOSUM62 similarity over the columns both sequences have, for rows of x."""
    num = np.zeros((len(x), len(x)))
    self_i = np.zeros_like(num)
    self_j = np.zeros_like(num)
    shared = np.zeros_like(num)
    for c in range(x.shape[1]):
        col = x[:, c]
        valid = col >= 0
        both = np.outer(valid, valid)
        safe = np.where(valid, col, 0)
        diag = B[safe, safe]
        num += np.where(both, B[safe][:, safe], 0)
        self_i += np.where(both, diag[:, None], 0)
        self_j += np.where(both, diag[None, :], 0)
        shared += both
    with np.errstate(invalid="ignore", divide="ignore"):
        d = 1 - num / np.sqrt(self_i * self_j)
    d[shared < args.min_shared * x.shape[1]] = np.nan
    np.fill_diagonal(d, 0)
    return d


def identity(x):
    """Share of identical residues over the columns both sequences have."""
    same = np.zeros((len(x), len(x)))
    shared = np.zeros_like(same)
    for c in range(x.shape[1]):
        col = x[:, c]
        valid = col >= 0
        both = np.outer(valid, valid)
        same += both & (col[:, None] == col[None, :])
        shared += both
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(shared > 0, same / shared, 0)


share = np.array([[bool(labels[i] & labels[j]) for j in range(n)] for i in range(n)])
whole = encode([aligned[sid] for sid in scored])
close_relative = identity(whole) >= args.max_identity
np.fill_diagonal(close_relative, True)
other_lineage = ~close_relative
has_other = other_lineage.any(axis=1)


def nearest_credit(d, allowed, shares):
    """Per transporter: share of its nearest allowed transporters that share a substrate."""
    d = np.where(allowed & ~np.isnan(d), d, np.inf)
    np.fill_diagonal(d, np.inf)
    best = d.min(axis=1, keepdims=True)
    nearest = np.isfinite(d) & (d <= best + 1e-9)
    counts = nearest.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        credit = np.where(counts > 0, (nearest & shares).sum(axis=1) / counts, 0.0)
    return credit, nearest


def knn_credit(d, k=3, seed=1):
    """(c): the substrates most common among the k nearest; credit for those that are its own."""
    rng = np.random.default_rng(seed)
    d = np.where(np.isnan(d), np.inf, d)
    credit = np.zeros(n)
    for i in range(n):
        order = sorted((j for j in range(n) if j != i), key=lambda j: (d[i, j], rng.random()))[:k]
        votes = {}
        for j in order:
            for s in labels[j]:
                votes[s] = votes.get(s, 0) + 1
        if votes:
            top = max(votes.values())
            predicted = [s for s, v in votes.items() if v == top]
            credit[i] = sum(s in labels[i] for s in predicted) / len(predicted)
    return credit


def shuffled_shares(d, allowed, rows):
    """Mean scores with labels shuffled over the transporters."""
    values, rng = [], random.Random(2)
    index = list(range(n))
    for _ in range(args.permutations):
        rng.shuffle(index)
        perm = share[np.ix_(index, index)]
        values.append(nearest_credit(d, allowed, perm)[0][rows].mean())
    return np.array(values)


def sign_flip_p(x, y, rows, flips=10000):
    diff = (x - y)[rows]
    observed = abs(diff.mean())
    rng = np.random.default_rng(3)
    signs = rng.choice([-1, 1], size=(flips, len(diff)))
    return (np.sum(np.abs((signs * diff).mean(axis=1)) >= observed - 1e-12) + 1) / (flips + 1)


# Distances of each method among the scored transporters
methods = {}
for shell in ("core", "extended"):
    methods[f"pocket_{shell}"] = blosum_distance(encode([signatures[sid][shell] for sid in scored]))
methods["all_columns"] = blosum_distance(whole[:, core_columns])
tree = Phylo.read(SEQUENCE_TREE, "newick")
tree.root_at_midpoint()
tips = {t.name: t for t in tree.get_terminals()}
lengths = path_lengths(tree, [tips[sid] for sid in scored])
methods["sequence_ml"] = np.array([[0.0 if a == b else lengths[(tips[a], tips[b])] for b in scored] for a in scored])

everyone = np.ones((n, n), dtype=bool)
rows_all = np.ones(n, dtype=bool)
results, per_transporter = {}, []
for method, d in methods.items():
    nn, nearest = nearest_credit(d, everyone, share)
    cross, cross_nearest = nearest_credit(d, other_lineage, share)
    knn = knn_credit(d)
    results[method] = dict(nn=nn, cross=cross, knn=knn,
                           nn_chance=shuffled_shares(d, everyone, rows_all),
                           cross_chance=shuffled_shares(d, other_lineage, has_other))
    for i, sid in enumerate(scored):
        per_transporter.append((method, sid, round(nn[i], 3), round(cross[i], 3) if has_other[i] else "",
                                round(knn[i], 3), ",".join(scored[j] for j in np.flatnonzero(nearest[i])),
                                ",".join(scored[j] for j in np.flatnonzero(cross_nearest[i]))))

# Random sets of well-aligned columns outside the pocket, of each pocket's size
candidates = [c for c in core_columns if c not in pocket_columns]
rng = random.Random(4)
random_rows, random_scores = [], {}
for shell in ("core", "extended"):
    size = len(signatures[scored[0]][shell])
    scores = []
    for k in range(args.random_sets):
        d = blosum_distance(whole[:, sorted(rng.sample(candidates, size))])
        nn_share = nearest_credit(d, everyone, share)[0].mean()
        cross_share = nearest_credit(d, other_lineage, share)[0][has_other].mean()
        scores.append((nn_share, cross_share))
        random_rows.append((shell, size, k + 1, round(nn_share, 4), round(cross_share, 4)))
    random_scores[f"pocket_{shell}"] = np.array(scores)

summary = []
for method, r in results.items():
    nn, cross = r["nn"].mean(), r["cross"][has_other].mean()
    row = [method, len(signatures[scored[0]][method.split("_")[1]]) if method.startswith("pocket") else
           (len(core_columns) if method == "all_columns" else ""), n, round(nn, 4),
           round(r["nn_chance"].mean(), 4), round((np.sum(r["nn_chance"] >= nn) + 1) / (args.permutations + 1), 4),
           int(has_other.sum()), round(cross, 4), round(r["cross_chance"].mean(), 4),
           round((np.sum(r["cross_chance"] >= cross) + 1) / (args.permutations + 1), 4), round(r["knn"].mean(), 4)]
    if method == "sequence_ml":
        row += ["", ""]
    else:
        row += [round(sign_flip_p(r["nn"], results["sequence_ml"]["nn"], rows_all), 4),
                round(sign_flip_p(r["cross"], results["sequence_ml"]["cross"], has_other), 4)]
    if method in random_scores:
        rs = random_scores[method]
        row += [round((np.sum(rs[:, 0] >= nn) + 1) / (len(rs) + 1), 4),
                round((np.sum(rs[:, 1] >= cross) + 1) / (len(rs) + 1), 4)]
    else:
        row += ["", ""]
    summary.append(row)

os.makedirs(EVAL_DIR, exist_ok=True)


def write(name, header, rows):
    with open(os.path.join(EVAL_DIR, name), "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


write("summary.tsv", ["method", "columns", "transporters", "nn_share", "nn_chance", "nn_p_shuffled",
                      "transporters_with_other_lineage", "cross_share", "cross_chance", "cross_p_shuffled",
                      "knn3_share", "nn_p_vs_sequence_ml", "cross_p_vs_sequence_ml", "nn_p_vs_random_columns",
                      "cross_p_vs_random_columns"], summary)
write("random_columns.tsv", ["pocket", "columns", "set", "nn_share", "cross_share"], random_rows)
write("per_transporter.tsv", ["method", "structure_id", "nn_credit", "cross_credit", "knn3_credit", "nearest",
                              "nearest_other_lineage"], per_transporter)

# Neighbour-joining trees of all models from the pocket signatures, for viewing and for clade purity
os.makedirs(TREE_DIR, exist_ok=True)
tree_paths = []
for shell in ("core", "extended"):
    d = blosum_distance(encode([signatures[sid][shell] for sid in models]))
    d = np.where(np.isnan(d), np.nanmax(d), d)
    matrix = [[float(d[i, j]) for j in range(i)] + [0.0] for i in range(len(models))]
    nj = DistanceTreeConstructor().nj(DistanceMatrix(models, matrix))
    for clade in nj.find_clades():
        if not clade.is_terminal():
            clade.name = None
        if clade.branch_length is not None and clade.branch_length < 0:
            clade.branch_length = 0.0
    nj.root_at_midpoint()
    path = os.path.join(TREE_DIR, f"pocket_{shell}.nwk")
    Phylo.write(nj, path, "newick")
    tree_paths.append(path)
subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "evaluate_substrate_clustering.py"), *tree_paths, SEQUENCE_TREE,
                "--out-dir", os.path.join(EVAL_DIR, "trees")], check=True, stdout=subprocess.DEVNULL)

for row in summary:
    print(f"  {row[0]}: (a) {row[3]:.1%} (chance {row[4]:.1%}); (a') other lineages {row[7]:.1%} of {row[6]} "
          f"(chance {row[8]:.1%}); (c) {row[10]:.1%}")
core = next(r for r in summary if r[0] == "pocket_core")
seq = next(r for r in summary if r[0] == "sequence_ml")
checks = [("(a) above the sequence tree", core[3] > seq[3] and core[11] < 0.05),
          ("(a') above the sequence tree", core[7] > seq[7] and core[12] < 0.05),
          ("(a) above 95% of random column sets", core[13] < 0.05),
          ("(a') above 95% of random column sets", core[14] < 0.05)]
for label, ok in checks:
    print(f"  core pocket {label}: {'yes' if ok else 'no'}")
print(f"Decision rule {'met' if all(ok for _, ok in checks) else 'not met'}. Results in {EVAL_DIR}/")
