#!/usr/bin/env python3
# clean_proteomes.py
# Description: Rewrites the collapsed proteomes with short protein IDs, <portal>-<proteinId>
#              (e.g. Dicsqu464_2-3083339), for alignments, trees, OrthoFinder and HMM searches.
#              MycoCosm headers ('jgi|portal|proteinId|name ...') contain '|', and some names
#              contain '#' or '_'-replaced '/', which break or get renamed in Newick trees.
#              Portal IDs only contain letters, digits and '_', so the last '-' always separates
#              the portal from the proteinId.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-20
# Version: 1.0
#
# Usage (from the project folder):
#   module load python-data
#   python3 scripts/clean_proteomes.py [portal_list] [proteomes_dir] [output_dir]
#   portal_list    TSV with a header row and portal IDs in column 1
#                  (default: 00_inputs/cohort_datasets.tsv).
#   proteomes_dir  Collapsed proteomes (default: 03_isoforms/collapsed_proteomes).
#   output_dir     Where the results are saved (default: 05_clean_proteomes/cohort).
# For the shadow proteomes:
#   python3 scripts/clean_proteomes.py 00_inputs/shadows_datasets.tsv \
#       03_isoforms/collapsed_proteomes 05_clean_proteomes/shadows
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

import os
import re
import sys

HEADER = re.compile(r"^>jgi\|([^|\s]+)\|(\d+)\|(\S*)")
MAP_COLUMNS = ["clean_id", "portal", "proteinId", "name", "original_header"]


def fail(portal, message):
    print(f"ERROR ({portal}): {message}", file=sys.stderr)
    sys.exit(1)


def read_portals(path):
    """Portal IDs from column 1 of a TSV with a header row."""
    with open(path) as f:
        lines = f.read().splitlines()[1:]
    portals = [line.split("\t")[0].strip() for line in lines]
    return [p for p in portals if p]


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


portal_list = sys.argv[1] if len(sys.argv) > 1 else "00_inputs/cohort_datasets.tsv"
proteomes_dir = sys.argv[2] if len(sys.argv) > 2 else "03_isoforms/collapsed_proteomes"
out_dir = sys.argv[3] if len(sys.argv) > 3 else "05_clean_proteomes/cohort"
usage = f"Usage: {sys.argv[0]} [portal_list] [proteomes_dir] [output_dir]"

if not os.path.isfile(portal_list):
    print(f"Portal list file '{portal_list}' not found\n{usage}", file=sys.stderr)
    sys.exit(1)
if not os.path.isdir(proteomes_dir):
    print(f"Proteomes folder '{proteomes_dir}' not found\n{usage}", file=sys.stderr)
    sys.exit(1)

portals = read_portals(portal_list)
missing = [p for p in portals if not os.path.isfile(os.path.join(proteomes_dir, f"{p}_proteins.faa"))]
if missing:
    print(f"{len(missing)} portals on the list have no proteome in '{proteomes_dir}':", file=sys.stderr)
    print("\n".join(f"  {p}" for p in missing), file=sys.stderr)
    sys.exit(1)

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
