#!/usr/bin/env python3
# map_pockets.py
# Description: Maps the sugar-binding pocket of three sugar-porter crystal structures onto every
#              sugar transporter in 08_structures/structure_ids.tsv, for the pocket concept test.
#              1. Pocket residues: residues of the reference with any heavy atom within 5 A (core)
#                 or 8 A (extended) of the bound sugar.
#              2. Mapping (a), through sequence: MAFFT L-INS-i alignment of all proteins with the
#                 three reference sequences; each pocket residue defines one alignment column.
#              3. Mapping (b), through structure: Foldseek alignment of each trimmed AlphaFold DB
#                 model to each reference structure.
#              A pocket column is kept when (a) and (b) give the same residue for >= min-agreement
#              of the models; the pocket signature of a protein is its residues in the kept columns.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-23
# Version: 1.0
#
# Usage (from the project folder; needs internet once, Biopython, mafft and foldseek):
#   python3 scripts/map_pockets.py [--min-agreement 0.8] [--threads N]
#
# References (sugar in the central binding site):
#   xyle   Escherichia coli XylE, PDB 4GBY (D-xylose) and 4GBZ (D-glucose), inward-facing occluded
#   glut3  Homo sapiens GLUT3, PDB 4ZW9 (D-glucose, both anomers), outward-occluded
#   stp10  Arabidopsis thaliana STP10, PDB 6H7D (D-glucose), outward-occluded
#
# Output, inside 08_structures/pockets:
#   references/<pdb>.cif, references/<name>.pdb   downloaded entries; protein chain of the first entry
#   reference_pockets.tsv       pocket residues of each reference, with their distance to the sugar
#   alignment.faa               MAFFT alignment of all proteins and the references (reused when
#                               its sequences are unchanged)
#   foldseek_references.tsv     Foldseek alignments of the models to the references
#   pocket_columns.tsv          every pocket column: shell, defining residues, (a)-(b) agreement, kept
#   mapping_agreement.tsv       per protein with a model: share of pocket columns where (a) and (b) agree
#   pocket_signatures.tsv       per protein: core and extended pocket residues in the kept columns
#   core_columns.txt            alignment columns with residues in >= 90% of the models (for baselines)

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import MMCIFParser, NeighborSearch, PDBIO, PDBParser, Select

IDS_TABLE = "08_structures/structure_ids.tsv"
MODELS_DIR = "08_structures/models_trimmed"
AFDB_DIR = "08_structures/afdb_models"
OUT_DIR = "08_structures/pockets"
REF_DIR = os.path.join(OUT_DIR, "references")
REFERENCES = [
    ("xyle", "Escherichia coli XylE", [("4GBY", {"XYP"}), ("4GBZ", {"BGC"})]),
    ("glut3", "Homo sapiens GLUT3", [("4ZW9", {"GLC", "BGC"})]),
    ("stp10", "Arabidopsis thaliana STP10", [("6H7D", {"BGC"})]),
]
SHELLS = [("core", 5.0), ("extended", 8.0)]

parser = argparse.ArgumentParser(description="Map sugar-binding pockets onto the sugar transporters")
parser.add_argument("--min-agreement", type=float, default=0.8)
parser.add_argument("--threads", type=int, default=os.cpu_count())
args = parser.parse_args()

for tool in ("mafft", "foldseek"):
    if shutil.which(tool) is None:
        sys.exit(f"ERROR: {tool} not found; activate the structure_test environment")
if not os.path.isdir(MODELS_DIR):
    sys.exit(f"ERROR: '{MODELS_DIR}' not found, run scripts/fetch_afdb_models.py first")
os.makedirs(REF_DIR, exist_ok=True)


def amino_acid(residue):
    """One-letter code of a protein residue (modified residues as their parent), or None."""
    if residue.id[0] not in (" ", "H_MSE") or "CA" not in residue:
        return None
    return protein_letters_3to1_extended.get(residue.get_resname())


def protein_residues(chain):
    return [r for r in chain if amino_acid(r)]


class Keep(Select):
    def __init__(self, residues):
        self.keep = {r.id for r in residues}

    def accept_residue(self, residue):
        return residue.id in self.keep


# 1. Reference pockets
references, pocket_rows = {}, []
for name, protein, entries in REFERENCES:
    distances = {}
    for n, (pdb_id, ligands) in enumerate(entries):
        path = os.path.join(REF_DIR, f"{pdb_id}.cif")
        if not os.path.isfile(path):
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.cif", path + ".part")
            os.replace(path + ".part", path)
        chain = MMCIFParser(QUIET=True).get_structure(pdb_id, path)[0]["A"]
        residues = protein_residues(chain)
        sugar = [a for r in chain if r.get_resname() in ligands for a in r if a.element != "H"]
        if not sugar:
            sys.exit(f"ERROR: no {'/'.join(sorted(ligands))} in {pdb_id} chain A")
        search = NeighborSearch([a for r in residues for a in r if a.element != "H"])
        for atom in sugar:
            for close in search.search(atom.coord, SHELLS[-1][1]):
                number = close.get_parent().id[1]
                d = float(close - atom)
                distances[number] = min(d, distances.get(number, d))
        if n == 0:
            # The first entry gives the reference sequence and the structure for Foldseek
            io = PDBIO()
            io.set_structure(chain)
            io.save(os.path.join(REF_DIR, f"{name}.pdb"), Keep(residues))
            numbers = [r.id[1] for r in residues]
            sequence = "".join(amino_acid(r) for r in residues)
    references[name] = dict(numbers=numbers, sequence=sequence, distances=distances)
    for number in sorted(distances):
        if number in numbers:
            shell = next(s for s, cut in SHELLS if distances[number] <= cut)
            pocket_rows.append((name, ",".join(e[0] for e in entries), number,
                                sequence[numbers.index(number)], round(distances[number], 2), shell))

# 2. Mapping (a): alignment of all proteins with the references
with open(IDS_TABLE) as f:
    proteins = [row for row in csv.DictReader(f, delimiter="\t") if row["use_in_structure_analysis"] == "TRUE"]
sequences = {p["structure_id"]: p["sequence"] for p in proteins}
alignment_path = os.path.join(OUT_DIR, "alignment.faa")
inputs = dict(list(sequences.items()) + [(f"ref_{n}", r["sequence"]) for n, r in references.items()])


def read_alignment(path):
    aligned, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                aligned[name] = []
            elif name:
                aligned[name].append(line.upper())
    return {k: "".join(v) for k, v in aligned.items()}


# An existing alignment of exactly these sequences is reused (L-INS-i takes minutes)
aligned = read_alignment(alignment_path) if os.path.isfile(alignment_path) else {}
if {k: v.replace("-", "") for k, v in aligned.items()} != inputs:
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as f:
        for sid, seq in inputs.items():
            f.write(f">{sid}\n{seq}\n")
        unaligned = f.name
    with open(alignment_path + ".part", "w") as out:
        subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", str(args.threads), "--quiet",
                        unaligned], stdout=out, check=True)
    os.replace(alignment_path + ".part", alignment_path)
    os.remove(unaligned)
    aligned = read_alignment(alignment_path)
n_columns = len(next(iter(aligned.values())))


def column_positions(gapped):
    """Alignment column -> 1-based position in the ungapped sequence (None at gaps)."""
    positions, position = [], 0
    for char in gapped:
        if char == "-":
            positions.append(None)
        else:
            position += 1
            positions.append(position)
    return positions


positions = {k: column_positions(v) for k, v in aligned.items()}

# Pocket columns: each pocket residue of a reference marks the column it is aligned in
columns = {}
for ref, refs, number, aa, distance, shell in pocket_rows:
    index = references[ref]["numbers"].index(number) + 1
    column = positions[f"ref_{ref}"].index(index)
    entry = columns.setdefault(column, {"shell": "extended", "defined_by": []})
    entry["defined_by"].append(f"{ref}:{aa}{number}")
    if shell == "core":
        entry["shell"] = "core"

# 3. Mapping (b): Foldseek alignments of each model to each reference
models = sorted(f[:-4] for f in os.listdir(MODELS_DIR) if f.endswith(".pdb"))
foldseek_path = os.path.join(OUT_DIR, "foldseek_references.tsv")
ref_structures = os.path.join(REF_DIR, "structures")
os.makedirs(ref_structures, exist_ok=True)
for name in references:
    shutil.copy(os.path.join(REF_DIR, f"{name}.pdb"), ref_structures)
with tempfile.TemporaryDirectory() as tmp:
    subprocess.run(["foldseek", "easy-search", MODELS_DIR, ref_structures, foldseek_path, tmp,
                    "--exhaustive-search", "1", "-e", "inf", "--threads", str(args.threads),
                    "--format-output", "query,target,qstart,tstart,qaln,taln,alntmscore"],
                   check=True, stdout=subprocess.DEVNULL)
shutil.rmtree(ref_structures)


def strip_name(name, known):
    for suffix in (".pdb_A", "_A", ".pdb"):
        if name.endswith(suffix) and name[:-len(suffix)] in known:
            return name[:-len(suffix)]
    return name


# Residue numbers of each model in file order; models are numbered as their AFDB sequence, which is
# our sequence for exact matches and is aligned to it for near-identical ones
aligner = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
model_numbers, model_sequences, to_ours = {}, {}, {}
for sid in models:
    chain = next(MMCIFParser(QUIET=True).get_structure(sid, os.path.join(AFDB_DIR, f"{sid}.cif"))[0].get_chains())
    afdb_sequence = "".join(amino_acid(r) or "X" for r in chain)
    afdb_numbers = [r.id[1] for r in chain]
    trimmed = next(PDBParser(QUIET=True).get_structure(sid, os.path.join(MODELS_DIR, f"{sid}.pdb"))[0].get_chains())
    model_numbers[sid] = [r.id[1] for r in protein_residues(trimmed)]
    model_sequences[sid] = "".join(amino_acid(r) for r in protein_residues(trimmed))
    if afdb_sequence == sequences[sid]:
        to_ours[sid] = {n: n for n in afdb_numbers}
    else:
        pair = aligner.align(afdb_sequence, sequences[sid])[0]
        to_ours[sid] = {}
        for (a0, a1), (b0, b1) in zip(*pair.aligned):
            for k in range(a1 - a0):
                to_ours[sid][afdb_numbers[a0 + k]] = b0 + k + 1

structural = {}  # (model, reference) -> {reference residue number: our position}
mismatches = 0
with open(foldseek_path) as f:
    for line in f:
        query, target, qstart, tstart, qaln, taln, _ = line.rstrip("\n").split("\t")
        query, target = strip_name(query, set(models)), strip_name(target, set(references))
        qi, ti, mapping = int(qstart) - 1, int(tstart) - 1, {}
        for q, t in zip(qaln, taln):
            if q != "-" and t != "-":
                mismatches += (q != model_sequences[query][qi]) + (t != references[target]["sequence"][ti])
                ours = to_ours[query].get(model_numbers[query][qi])
                if ours:
                    mapping[references[target]["numbers"][ti]] = ours
            qi += q != "-"
            ti += t != "-"
        structural[(query, target)] = mapping
if mismatches:
    sys.exit(f"ERROR: {mismatches} aligned residues in Foldseek's output do not match the structures' residues")
missing = [(m, r) for m in models for r in references if (m, r) not in structural]
if missing:
    sys.exit(f"ERROR: Foldseek gave no alignment for {len(missing)} model-reference pairs, e.g. {missing[0]}")

# Agreement per column: (a) and (b) give the same residue of a model for at least one reference
# that has a residue in that column (a reference residue in the column, through (b), points to the
# residue the alignment puts there)
agree = {}
for column in columns:
    for sid in models:
        ours = positions[sid][column]
        hits = []
        for ref in references:
            index = positions[f"ref_{ref}"][column]
            if index is None:
                continue
            hits.append(structural[(sid, ref)].get(references[ref]["numbers"][index - 1]) == ours)
        agree[(column, sid)] = ours is not None and any(hits)
column_rows, kept = [], {"core": [], "extended": []}
for column in sorted(columns):
    with_residue = [sid for sid in models if positions[sid][column] is not None]
    rate = sum(agree[(column, sid)] for sid in with_residue) / len(with_residue) if with_residue else 0.0
    keep = rate >= args.min_agreement
    shell = columns[column]["shell"]
    if keep:
        kept["extended"].append(column)
        if shell == "core":
            kept["core"].append(column)
    column_rows.append((column + 1, shell, ";".join(columns[column]["defined_by"]), len(with_residue),
                        round(rate, 3), keep))

agreement_rows = []
for sid in models:
    row = [sid]
    for shell in ("core", "extended"):
        shell_columns = [c for c in columns if shell == "extended" or columns[c]["shell"] == "core"]
        row.append(round(sum(agree[(c, sid)] for c in shell_columns) / len(shell_columns), 3))
    agreement_rows.append(row)

core_columns = [c for c in range(n_columns) if sum(aligned[sid][c] != "-" for sid in models) >= 0.9 * len(models)]


def write(name, header, rows):
    with open(os.path.join(OUT_DIR, name), "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


write("reference_pockets.tsv", ["reference", "pdb_entries", "residue", "amino_acid", "min_distance_A", "shell"],
      pocket_rows)
write("pocket_columns.tsv", ["column", "shell", "defined_by", "models_with_residue", "agreement", "kept"],
      column_rows)
write("mapping_agreement.tsv", ["structure_id", "core_agreement", "extended_agreement"], agreement_rows)
write("pocket_signatures.tsv", ["structure_id", "has_model", "core", "extended"],
      [(sid, sid in model_numbers, "".join(aligned[sid][c] for c in kept["core"]),
        "".join(aligned[sid][c] for c in kept["extended"])) for sid in sequences])
with open(os.path.join(OUT_DIR, "core_columns.txt"), "w") as f:
    f.write("\n".join(str(c + 1) for c in core_columns) + "\n")

n_core = sum(1 for c in columns if columns[c]["shell"] == "core")
print(f"Pocket columns: {n_core} core and {len(columns) - n_core} extended-only; kept after (a)-(b) agreement "
      f">= {args.min_agreement:.0%}: {len(kept['core'])} core, {len(kept['extended'])} extended. "
      f"{len(core_columns)} alignment columns have residues in >= 90% of the {len(models)} models. "
      f"Results in {OUT_DIR}/")
