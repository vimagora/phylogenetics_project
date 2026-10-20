#!/usr/bin/env python3
# collect_marker_trees.py
# Description: Collects the per-marker trees and best-fit models (marker_trees.sbatch): writes all
#              trees into one file for ASTRAL, and a partition file with each marker's model for
#              the concatenated IQ-TREE run, so it does not repeat the model selection.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-20
# Version: 1.0
#
# Usage (from the project folder):
#   module load python-data
#   python3 scripts/collect_marker_trees.py [trees_dir] [concatenated_dir] [output_dir]
#   trees_dir         Per-marker IQ-TREE results (default: 06_species_tree/marker_trees).
#   concatenated_dir  Output of concatenate_markers.py, with kept_markers.txt and partitions.nex
#                     (default: 06_species_tree/concatenated).
#   output_dir        Where the results are saved (default: 06_species_tree).
#
# Output, inside output_dir:
#   marker_trees.nwk             one tree per line, in the order of kept_markers.txt
#   marker_models.tsv            best-fit model (BIC) of each marker
#   concatenated/partitions_models.nex  partitions.nex plus a charpartition with each marker's model
#                                       (written next to partitions.nex)
# The run stops if a kept marker has no tree or no best-fit model.

import os
import re
import sys

trees_dir = sys.argv[1] if len(sys.argv) > 1 else "06_species_tree/marker_trees"
concat_dir = sys.argv[2] if len(sys.argv) > 2 else "06_species_tree/concatenated"
out_dir = sys.argv[3] if len(sys.argv) > 3 else "06_species_tree"
BEST_MODEL = re.compile(r"^Best-fit model according to BIC:\s+(\S+)", re.MULTILINE)
CHARSET = re.compile(r"^\s*charset\s+(\S+)\s*=\s*(\d+)-(\d+);", re.MULTILINE)


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


kept_path = os.path.join(concat_dir, "kept_markers.txt")
partitions_path = os.path.join(concat_dir, "partitions.nex")
for path in (trees_dir, kept_path, partitions_path):
    if not os.path.exists(path):
        fail(f"'{path}' not found")

with open(kept_path) as f:
    markers = [line.strip() for line in f if line.strip()]
with open(partitions_path) as f:
    charsets = CHARSET.findall(f.read())
if [name for name, _, _ in charsets] != markers:
    fail(f"the charsets in '{partitions_path}' do not match '{kept_path}'")

trees, models, missing = [], [], []
for marker in markers:
    tree_path = os.path.join(trees_dir, f"{marker}.treefile")
    report_path = os.path.join(trees_dir, f"{marker}.iqtree")
    if not (os.path.isfile(tree_path) and os.path.isfile(report_path)):
        missing.append(marker)
        continue
    with open(tree_path) as f:
        trees.append(f.read().strip())
    with open(report_path) as f:
        match = BEST_MODEL.search(f.read())
    if not match:
        fail(f"no best-fit model in '{report_path}'")
    models.append(match.group(1))
if missing:
    fail(f"{len(missing)} markers have no finished tree, e.g. {missing[0]}; "
         f"resubmit scripts/marker_trees.sbatch")

with open(os.path.join(out_dir, "marker_trees.nwk"), "w") as f:
    f.write("\n".join(trees) + "\n")
with open(os.path.join(out_dir, "marker_models.tsv"), "w") as f:
    f.write("busco_id\tbest_model_bic\n")
    f.write("".join(f"{m}\t{model}\n" for m, model in zip(markers, models)))
with open(os.path.join(concat_dir, "partitions_models.nex"), "w") as f:
    f.write("#nexus\nbegin sets;\n")
    for name, start, end in charsets:
        f.write(f"    charset {name} = {start}-{end};\n")
    f.write("    charpartition markers = " + ", ".join(f"{model}:{m}" for m, model in zip(markers, models)) + ";\n")
    f.write("end;\n")

counts = {}
for model in models:
    base = model.split("+")[0]
    counts[base] = counts.get(base, 0) + 1
print(f"Collected {len(trees)} marker trees into {out_dir}/marker_trees.nwk")
print("Best-fit base models: " + ", ".join(f"{b} {n}" for b, n in sorted(counts.items(), key=lambda x: -x[1])))
print(f"Partition file with models: {concat_dir}/partitions_models.nex")
