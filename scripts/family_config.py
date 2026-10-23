#!/usr/bin/env python3
# family_config.py
# Description: Shared helpers of the generic family scripts: reads a family's settings
#              (families/<family>/family.conf) and the tables every family script needs.
#              Not run on its own; imported by summarize_family.py, build_family_fasta.py,
#              fetch_outgroups.py and make_family_itol.py.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-23
# Version: 1.0

import csv
import os
import re
import sys

TAXONOMY = "00_inputs/cohort_taxonomy.tsv"
OG_DIR = "07_orthofinder/run/Results_cohort/Orthogroups"


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
    """Cohort taxonomy by portal (Windows line endings tolerated)."""
    return {r["portal"].strip(): {k: v.strip() for k, v in r.items()} for r in read_tsv(TAXONOMY)}


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
