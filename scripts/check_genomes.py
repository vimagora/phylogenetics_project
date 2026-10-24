#!/usr/bin/env python3
# check_genomes.py
# Description: Checks 00_inputs/genomes.tsv before the second half of the pipeline: unique portals,
#              known roles (cohort, shadow), taxonomy and NCBI taxon ID filled for every cohort
#              genome, every genome listed in candidates.tsv, and a collapsed proteome and BUSCO
#              results for every genome. clean_proteomes.py runs the same checks; this script lets
#              you run them on their own while choosing the cohort.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-24
# Version: 1.0
#
# Usage (from the project folder): python3 scripts/check_genomes.py [--genomes FILE] [--role cohort|shadow]

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genome_tables  # noqa: E402

parser = argparse.ArgumentParser(description="Check genomes.tsv before the second half of the pipeline")
parser.add_argument("--genomes", default=genome_tables.GENOMES)
parser.add_argument("--role", choices=genome_tables.ROLES, help="check only this role (default: all)")
args = parser.parse_args()

roles = (args.role,) if args.role else genome_tables.ROLES
genome_tables.require(args.genomes, roles)
counts = {r: len(genome_tables.portals(args.genomes, (r,))) for r in roles}
print(f"{args.genomes} is ready: " + ", ".join(f"{n} {r}" for r, n in counts.items()))
