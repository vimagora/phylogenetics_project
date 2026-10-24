#!/usr/bin/env python3
# clean_proteomes.py
# Description: Rewrites the collapsed proteomes with short protein IDs, <portal>-<proteinId>
#              (e.g. Dicsqu464_2-3083339), for alignments, trees, OrthoFinder and HMM searches.
#              MycoCosm headers ('jgi|portal|proteinId|name ...') contain '|', and some names
#              contain '#' or '_'-replaced '/', which break or get renamed in Newick trees.
#              Portal IDs only contain letters, digits and '_', so the last '-' always separates
#              the portal from the proteinId.
# Author: Victor M. Gonzalez Ramos
#              It is the first step of the second half of the pipeline, so it first checks
#              00_inputs/genomes.tsv (see check_genomes.py) and stops if a genome is not ready.
# Date: 2026-10-20
# Version: 1.1
#
# Usage (from the project folder):
#   module load python-data
#   python3 scripts/clean_proteomes.py [--role cohort|shadow] [--genomes FILE] [--proteomes-dir DIR]
#       [--out-dir DIR]
#   --role           Genomes of this role in genomes.tsv (default: cohort).
#   --proteomes-dir  Collapsed proteomes (default: 03_isoforms/collapsed_proteomes).
#   --out-dir        Where the results are saved (default: 05_clean_proteomes/cohort or .../shadows).
#
# Output, inside output_dir:
#   proteomes/<portal>.faa  one protein per line, header '>portal-proteinId' only; the file name
#                           is the portal, which OrthoFinder uses as the species name
#   id_map.tsv              clean_id, portal, proteinId, gene name and the original header
#
# Changes to the sequences: the terminal '*' (stop codon) is removed, and any internal '*' is
# replaced by 'X' and counted. The run stops if a header is not 'jgi|<portal>|<proteinId>|...'
# with the portal of the file name, if a proteinId appears twice, if a sequence is empty, or if a
# sequence has characters other than letters. All proteomes are rewritten on every run.

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genome_tables  # noqa: E402

HEADER = re.compile(r"^>jgi\|([^|\s]+)\|(\d+)\|(\S*)")
MAP_COLUMNS = ["clean_id", "portal", "proteinId", "name", "original_header"]


def fail(portal, message):
    print(f"ERROR ({portal}): {message}", file=sys.stderr)
    sys.exit(1)


def read_fasta(path):
    """(header, sequence) pairs; ';' comment lines are skipped."""
    records, header, chunks = [], None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line, []
            else:
                chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def clean_proteome(portal, in_path, out_path):
    """Writes the cleaned proteome and returns its id_map rows and the number of internal '*'."""
    rows, seen, n_internal = [], set(), 0
    with open(out_path + ".part", "w") as out:
        for header, seq in read_fasta(in_path):
            match = HEADER.match(header)
            if not match:
                fail(portal, f"header not in 'jgi|portal|proteinId|name' format: {header}")
            header_portal, protein_id, name = match.groups()
            if header_portal != portal:
                fail(portal, f"header of another portal ({header_portal}): {header}")
            if protein_id in seen:
                fail(portal, f"proteinId {protein_id} appears twice")
            seen.add(protein_id)
            seq = seq.upper().removesuffix("*")
            n_internal += seq.count("*")
            seq = seq.replace("*", "X")
            if not seq:
                fail(portal, f"empty sequence: {header}")
            if not seq.isalpha():
                fail(portal, f"unexpected characters in sequence: {header}")
            clean_id = f"{portal}-{protein_id}"
            out.write(f">{clean_id}\n{seq}\n")
            rows.append((clean_id, portal, protein_id, name, header[1:].replace("\t", " ")))
    if not rows:
        fail(portal, f"no sequences in '{in_path}'")
    os.replace(out_path + ".part", out_path)
    return rows, n_internal


parser = argparse.ArgumentParser(description="Clean protein IDs of the collapsed proteomes")
parser.add_argument("--role", choices=genome_tables.ROLES, default="cohort")
parser.add_argument("--genomes", default=genome_tables.GENOMES)
parser.add_argument("--proteomes-dir", default=genome_tables.COLLAPSED_DIR)
parser.add_argument("--out-dir")
args = parser.parse_args()
proteomes_dir = args.proteomes_dir
out_dir = args.out_dir or os.path.join("05_clean_proteomes", {"cohort": "cohort", "shadow": "shadows"}[args.role])

genome_tables.require(args.genomes, (args.role,), collapsed_dir=proteomes_dir)
portals = genome_tables.portals(args.genomes, (args.role,))

os.makedirs(os.path.join(out_dir, "proteomes"), exist_ok=True)
map_path = os.path.join(out_dir, "id_map.tsv")
n_proteins = 0
with open(map_path + ".part", "w") as id_map:
    id_map.write("\t".join(MAP_COLUMNS) + "\n")
    for portal in portals:
        rows, n_internal = clean_proteome(portal, os.path.join(proteomes_dir, f"{portal}_proteins.faa"),
                                          os.path.join(out_dir, "proteomes", f"{portal}.faa"))
        for row in rows:
            id_map.write("\t".join(row) + "\n")
        n_proteins += len(rows)
        note = f", {n_internal} internal '*' replaced by 'X'" if n_internal else ""
        print(f"  {portal}: {len(rows)} proteins{note}")
os.replace(map_path + ".part", map_path)

print(f"Cleaned {len(portals)} proteomes ({n_proteins} proteins) into {out_dir}/proteomes, "
      f"ID map: {map_path}")
