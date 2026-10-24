#!/usr/bin/env python3
# genome_tables.py
# Description: Shared helpers for the two genome tables (not run on its own):
#              00_inputs/candidates.tsv  every genome assessed in the first half of the pipeline
#                                        (download, assembly statistics, isoform collapse, BUSCO):
#                                        portal, ncbi_taxid, mycocosm_name, notes;
#              00_inputs/genomes.tsv     the genomes chosen for the analysis, with their role:
#                                        portal, ncbi_taxid, species, phylum, subphylum, class,
#                                        order, role (cohort or shadow), notes.
#              check() validates genomes.tsv and confirms that every chosen genome was processed
#              by the first half, so that a mistake stops the second half within seconds.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-24
# Version: 1.0

import csv
import os
import sys

CANDIDATES = "00_inputs/candidates.tsv"
GENOMES = "00_inputs/genomes.tsv"
COLLAPSED_DIR = "03_isoforms/collapsed_proteomes"
BUSCO_TABLE = "04_BUSCO/busco_full_table_collapsed_proteomes.tsv"
ROLES = ("cohort", "shadow")
RANKS = ("phylum", "subphylum", "class", "order")
GENOME_COLUMNS = ("portal", "ncbi_taxid", "species", *RANKS, "role", "notes")


def _fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def _read(path):
    if not os.path.isfile(path):
        _fail(f"'{path}' not found")
    with open(path, newline="", encoding="utf-8") as f:
        return [{k.strip(): (v or "").strip() for k, v in row.items()} for row in csv.DictReader(f, delimiter="\t")]


def candidates(path=CANDIDATES):
    """Candidate portals, in table order."""
    return [r["portal"] for r in _read(path) if r["portal"]]


def genomes(path=GENOMES, roles=("cohort",)):
    """Rows of genomes.tsv with the given roles, in table order."""
    return [r for r in _read(path) if r["portal"] and r["role"] in roles]


def portals(path=GENOMES, roles=("cohort",)):
    return [r["portal"] for r in genomes(path, roles)]


def check(path=GENOMES, roles=ROLES, collapsed_dir=COLLAPSED_DIR, busco_table=BUSCO_TABLE,
          candidates_path=CANDIDATES):
    """Problems in genomes.tsv for the given roles, as a list of messages (empty if none)."""
    rows = _read(path)
    problems = []
    missing_columns = [c for c in GENOME_COLUMNS if rows and c not in rows[0]]
    if missing_columns:
        return [f"{path} lacks the columns {', '.join(missing_columns)}"]
    seen = set()
    for r in rows:
        if not r["portal"]:
            continue
        if r["portal"] in seen:
            problems.append(f"{r['portal']} appears twice")
        seen.add(r["portal"])
        if r["role"] not in ROLES:
            problems.append(f"{r['portal']}: role '{r['role']}' is not one of {', '.join(ROLES)}")
        if r["role"] == "cohort":
            empty = [c for c in ("ncbi_taxid", "species", *RANKS) if not r[c]]
            if empty:
                problems.append(f"{r['portal']} (cohort): empty {', '.join(empty)}")
    chosen = [r["portal"] for r in rows if r["portal"] and r["role"] in roles]
    if not chosen:
        problems.append(f"no genomes with role {' or '.join(roles)} in {path}")
    known = set(candidates(candidates_path))
    problems += [f"{p} is not in {candidates_path}" for p in chosen if p not in known]
    problems += [f"{p} has no collapsed proteome in {collapsed_dir}" for p in chosen
                 if not os.path.isfile(os.path.join(collapsed_dir, f"{p}_proteins.faa"))]
    if not os.path.isfile(busco_table):
        problems.append(f"{busco_table} not found (scripts/combine_busco_tables.sh)")
    else:
        with open(busco_table) as f:
            in_busco = {line.split("\t", 1)[0] for line in f}
        problems += [f"{p} has no BUSCO results in {busco_table}" for p in chosen if p not in in_busco]
    return problems


def require(path=GENOMES, roles=ROLES, **paths):
    """Stops with a list of problems if check() finds any."""
    problems = check(path, roles, **paths)
    if problems:
        print(f"ERROR: {path} is not ready ({len(problems)} problems):", file=sys.stderr)
        print("\n".join(f"  {p}" for p in problems[:30]), file=sys.stderr)
        if len(problems) > 30:
            print(f"  ... and {len(problems) - 30} more", file=sys.stderr)
        sys.exit(1)
