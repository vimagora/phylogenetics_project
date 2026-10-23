#!/usr/bin/env python3
# fetch_afdb_models.py
# Description: Finds an AlphaFold Database (AFDB) model for each sugar transporter in
#              08_structures/structure_ids.tsv, downloads it, and trims its low-confidence termini.
#              A model is accepted when its sequence is identical to ours (found through UniParc by
#              sequence checksum), or, failing that, when a UniProt entry of the same species has
#              >= MIN_IDENTITY identity over the full length (e.g. another strain or an older version
#              of the same protein) and has an AFDB model.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-22
# Version: 1.0
#
# Usage (from the project folder; needs internet and Biopython):
#   python3 scripts/fetch_afdb_models.py [--min-identity 0.98] [--min-plddt 70]
#
# Output, inside 08_structures:
#   afdb_models/<structure_id>.cif      AFDB model as downloaded (pLDDT in the B-factor column)
#   models_trimmed/<structure_id>.pdb   the same model without its low-confidence termini: residues
#                                       before the first and after the last residue with pLDDT >=
#                                       min-plddt are removed; internal loops are kept
#   models.tsv                          per protein: match type, AFDB accession, identity, mean
#                                       pLDDT and the residues kept
# Downloaded models are reused when the script is run again.

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, PDBIO, Select

IDS_TABLE = "08_structures/structure_ids.tsv"
MODELS_DIR = "08_structures/afdb_models"
TRIMMED_DIR = "08_structures/models_trimmed"
SUMMARY = "08_structures/models.tsv"
UNIPROT = "https://rest.uniprot.org"
AFDB = "https://alphafold.ebi.ac.uk/api/prediction"


def get(url, retries=3):
    """Body of a URL as text, with a few retries; None on HTTP 404."""
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return response.read().decode()
        except urllib.error.HTTPError as error:
            if error.code in (400, 404):
                return None
            time.sleep(5 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError):
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"could not fetch {url}")


def afdb_entry(accession):
    """AFDB prediction record for a UniProt accession, or None."""
    body = get(f"{AFDB}/{accession}")
    records = json.loads(body) if body else []
    return records[0] if records else None


def exact_match(sequence):
    """AFDB record whose sequence is identical to ours, through UniParc's MD5 checksum."""
    md5 = hashlib.md5(sequence.encode()).hexdigest().upper()
    query = urllib.parse.urlencode({"query": f"checksum:{md5}", "fields": "upi,accession", "format": "tsv"})
    body = get(f"{UNIPROT}/uniparc/search?{query}") or ""
    accessions = [a.strip().split(".")[0] for line in body.splitlines()[1:]
                  for a in line.split("\t")[1].split(";") if a.strip()]
    for accession in accessions:
        entry = afdb_entry(accession)
        if entry and entry.get("sequence") == sequence:
            return entry
    return None


def species_taxon(species):
    """NCBI taxon ID of a species name such as 'Saccharomyces cerevisiae (synonym ...)'."""
    name = re.sub(r"\s*\(.*\)", "", species.split(" / ")[0]).strip()
    name = " ".join(name.split()[:2])
    query = urllib.parse.urlencode({"query": f'"{name}"', "fields": "id,scientific_name,rank", "format": "tsv",
                                    "size": 100})
    body = get(f"{UNIPROT}/taxonomy/search?{query}") or ""
    rows = [line.split("\t") for line in body.splitlines()[1:] if line.count("\t") >= 2]
    # The species itself (a plain text search also returns e.g. viruses named after it), else a strain
    for taxon, scientific_name, rank in rows:
        if rank == "species" and scientific_name.lower() == name.lower():
            return taxon
    for taxon, scientific_name, rank in rows:
        if rank == "strain" and scientific_name.lower().startswith(name.lower()):
            return taxon
    return None


def identity(a, b, aligner):
    """Identical positions in a global alignment, divided by the longer sequence's length."""
    alignment = aligner.align(a, b)[0]
    matches = sum(x == y != "-" for x, y in zip(alignment[0], alignment[1]))
    return matches / max(len(a), len(b))


def kmers(sequence, k=5):
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1)}


def near_identical_match(sequence, species, min_identity, aligner):
    """AFDB record of the most similar UniProt entry of the same species, if identity >= min_identity."""
    taxon = species_taxon(species)
    if taxon is None:
        return None, 0.0
    low, high = int(len(sequence) * 0.95), int(len(sequence) * 1.05) + 1
    query = urllib.parse.urlencode({"query": f"taxonomy_id:{taxon} AND length:[{low} TO {high}]", "format": "fasta"})
    body = get(f"{UNIPROT}/uniprotkb/stream?{query}") or ""
    candidates, accession = {}, None
    for line in body.splitlines():
        if line.startswith(">"):
            accession = line.split("|")[1]
            candidates[accession] = []
        elif accession:
            candidates[accession].append(line.strip())
    own = kmers(sequence)
    ranked = sorted(((len(own & kmers("".join(s))) / len(own), a, "".join(s)) for a, s in candidates.items()),
                    reverse=True)[:5]
    best = (None, 0.0)
    for shared, accession, candidate in ranked:
        if shared < 0.5:
            break
        score = identity(sequence, candidate, aligner)
        if score >= min_identity and score > best[1]:
            entry = afdb_entry(accession)
            if entry and entry.get("sequence") == candidate:
                best = (entry, score)
    return best


class KeepRange(Select):
    def __init__(self, keep):
        self.keep = keep

    def accept_residue(self, residue):
        return residue.id in self.keep


def trim(cif_path, pdb_path, min_plddt):
    """Removes low-confidence termini; returns (mean pLDDT of kept residues, first, last, n kept)."""
    structure = MMCIFParser(QUIET=True).get_structure("model", cif_path)
    residues = [r for r in structure[0].get_residues() if "CA" in r]
    plddt = [r["CA"].get_bfactor() for r in residues]
    confident = [i for i, value in enumerate(plddt) if value >= min_plddt]
    if not confident:
        return None
    first, last = confident[0], confident[-1]
    keep = {residues[i].id for i in range(first, last + 1)}
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path, KeepRange(keep))
    kept = plddt[first:last + 1]
    return round(sum(kept) / len(kept), 1), first + 1, last + 1, len(kept)


parser = argparse.ArgumentParser(description="Download and trim AlphaFold Database models")
parser.add_argument("--min-identity", type=float, default=0.98)
parser.add_argument("--min-plddt", type=float, default=70)
args = parser.parse_args()

if not os.path.isfile(IDS_TABLE):
    print(f"ERROR: '{IDS_TABLE}' not found, run scripts/make_af3_inputs.py first", file=sys.stderr)
    sys.exit(1)

aligner = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
os.makedirs(MODELS_DIR, exist_ok=True)
os.makedirs(TRIMMED_DIR, exist_ok=True)

with open(IDS_TABLE) as f:
    proteins = [row for row in csv.DictReader(f, delimiter="\t") if row["use_in_structure_analysis"] == "TRUE"]

rows = []
for protein in proteins:
    sid, sequence = protein["structure_id"], protein["sequence"]
    cif_path = os.path.join(MODELS_DIR, f"{sid}.cif")
    entry, match, score = exact_match(sequence), "exact", 1.0
    if entry is None:
        entry, score = near_identical_match(sequence, protein["species"], args.min_identity, aligner)
        match = "near_identical" if entry else "none"
    row = dict(structure_id=sid, match=match, afdb_accession="", identity="", query_length=len(sequence),
               model_length="", mean_plddt_trimmed="", first_kept="", last_kept="", residues_kept="")
    if entry:
        if not os.path.isfile(cif_path):
            with open(cif_path + ".part", "w") as f:
                f.write(get(entry["cifUrl"]))
            os.replace(cif_path + ".part", cif_path)
        trimmed = trim(cif_path, os.path.join(TRIMMED_DIR, f"{sid}.pdb"), args.min_plddt)
        row.update(afdb_accession=entry["uniprotAccession"], identity=round(score, 3),
                   model_length=len(entry["sequence"]))
        if trimmed:
            row.update(mean_plddt_trimmed=trimmed[0], first_kept=trimmed[1], last_kept=trimmed[2],
                       residues_kept=trimmed[3])
    rows.append(row)
    print(f"  {sid}: {match}" + (f" ({row['afdb_accession']}, identity {row['identity']})" if entry else ""))

with open(SUMMARY, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

counts = {m: sum(r["match"] == m for r in rows) for m in ("exact", "near_identical", "none")}
print(f"{len(rows)} proteins: {counts['exact']} exact AFDB models, {counts['near_identical']} near-identical "
      f"(>= {args.min_identity:.0%}), {counts['none']} without a model. Trimmed models in {TRIMMED_DIR}, "
      f"summary in {SUMMARY}")
