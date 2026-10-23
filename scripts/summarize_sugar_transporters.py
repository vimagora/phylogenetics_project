#!/usr/bin/env python3
# summarize_sugar_transporters.py
# Description: Turns the searches of search_sugar_transporters.sbatch into the cohort sugar
#              transporters and their orthogroups:
#              1. anchors: the cohort copy of each reference transporter, i.e. a DIAMOND hit with
#                 >= min-identity identity covering >= min-coverage of both proteins, in the
#                 reference's own genome (its listed portal, or the cohort genome of its species).
#                 A listed clean ID (D. squalens, some A. niger rows) needs only the query coverage,
#                 since the reference may be a trimmed gene model (DsST24). References without an
#                 anchor keep their best hit and its orthogroup.
#              2. sugar porters: cohort proteins with PF00083 (Sugar_tr) at the gathering threshold;
#                 those whose domains cover < min-hmm-coverage of the HMM are flagged as partial.
#              3. orthogroups (OrthoFinder Orthogroups.tsv, i.e. root hierarchical orthogroups):
#                 size, species, sugar porters and anchors of each orthogroup holding a sugar porter
#                 or an anchor.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-23
# Version: 1.0
#
# Usage (from the project folder; standard library only):
#   python3 scripts/summarize_sugar_transporters.py [--min-identity 98] [--min-coverage 0.9]
#       [--min-hmm-coverage 0.5] [--orthofinder-dir DIR]
#   --orthofinder-dir  OrthoFinder results (default: 07_orthofinder/run/Results_cohort).
#
# Output, inside 09_sugar_transporters:
#   anchors.tsv               per reference: status, anchor or best hit, identity, coverages, orthogroup
#   pf00083_hits.tsv          per cohort protein with PF00083: length, HMM coverage, score, partial, orthogroup
#   pf00083_per_genome.tsv    per genome: full-length and partial sugar porters
#   orthogroup_summary.tsv    per orthogroup with a sugar porter or an anchor
#   st_proteins.tsv           full-length sugar porters with orthogroup and anchor status (tree input)

import argparse
import csv
import os
import re
import sys

IDS_TABLE = "08_structures/structure_ids.tsv"
TAXONOMY = "00_inputs/cohort_taxonomy.tsv"
SEARCH_DIR = "09_sugar_transporters/search"
OUT_DIR = "09_sugar_transporters"
# Older names in the reference table for species named differently in the cohort taxonomy
SYNONYMS = {"piriformospora indica": "serendipita indica", "postia placenta": "rhodonia placenta",
            "myceliophthora thermophila": "thermothelomyces thermophilus",
            "pichia stipitis": "scheffersomyces stipitis", "magnaporthe oryzae": "pyricularia oryzae",
            "dekkera bruxellensis": "brettanomyces bruxellensis"}

parser = argparse.ArgumentParser(description="Cohort sugar transporters, anchors and orthogroups")
parser.add_argument("--min-identity", type=float, default=98)
parser.add_argument("--min-coverage", type=float, default=0.9)
parser.add_argument("--min-hmm-coverage", type=float, default=0.5)
parser.add_argument("--orthofinder-dir", default="07_orthofinder/run/Results_cohort")
args = parser.parse_args()

inputs = {name: os.path.join(SEARCH_DIR, name) for name in
          ("diamond_references.tsv", "pf00083_cohort.domtbl", "pf00083_references.domtbl")}
og_table = os.path.join(args.orthofinder_dir, "Orthogroups", "Orthogroups.tsv")
unassigned_table = os.path.join(args.orthofinder_dir, "Orthogroups", "Orthogroups_UnassignedGenes.tsv")
for path in list(inputs.values()) + [og_table, unassigned_table, IDS_TABLE, TAXONOMY]:
    if not os.path.isfile(path):
        sys.exit(f"ERROR: '{path}' not found (searches: scripts/search_sugar_transporters.sbatch)")


def binomial(name):
    return " ".join(re.sub(r"[^A-Za-z. -]", " ", name).split()[:2]).lower()


with open(TAXONOMY) as f:
    taxonomy = {r["portal"]: r for r in csv.DictReader(f, delimiter="\t")}
portals_of_species = {}
for portal, row in taxonomy.items():
    portals_of_species.setdefault(binomial(row["species"]), set()).add(portal)

with open(IDS_TABLE) as f:
    references = [r for r in csv.DictReader(f, delimiter="\t") if r["use_in_structure_analysis"] == "TRUE"]


def expected_portals(ref):
    """The reference's own cohort genome: its listed portal, else the cohort genome(s) of its species."""
    if ref.get("cohort_portal"):
        return {ref["cohort_portal"]}
    names = set()
    for part in ref["species"].split(" / "):
        names.add(binomial(part.split("(")[0]))
        for alias in re.findall(r"\((?:synonym|anamorph)\s+([^)]+)\)", part):
            names.add(binomial(alias))
    names |= {SYNONYMS[n] for n in names if n in SYNONYMS}
    return set().union(*(portals_of_species.get(n, set()) for n in names))


# Orthogroup of every cohort gene; unassigned genes are single-gene orthogroups
orthogroup, members = {}, {}
for path in (og_table, unassigned_table):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            genes = [g for cell in fields[1:] for g in cell.split(", ") if g]
            members[fields[0]] = genes
            for gene in genes:
                orthogroup[gene] = fields[0]


def portal_of(clean_id):
    return clean_id.rsplit("-", 1)[0]


# 1. Anchors
hits = {}
with open(inputs["diamond_references.tsv"]) as f:
    for line in f:
        q, s, pident, length, qlen, slen, qs, qe, ss, se, evalue, bits = line.rstrip("\n").split("\t")
        hits.setdefault(q, []).append(dict(subject=s, identity=float(pident), qcov=(int(qe) - int(qs) + 1) / int(qlen),
                                           scov=(int(se) - int(ss) + 1) / int(slen), bits=float(bits)))
anchor_rows, anchor_of = [], {}
for ref in references:
    sid, own = ref["structure_id"], expected_portals(ref)
    candidates = sorted(hits.get(sid, []), key=lambda h: -h["bits"])
    listed = ref.get("cohort_id", "")
    chosen, status = None, "no_hit"
    if listed:
        match = [h for h in candidates if h["subject"] == listed]
        if match and match[0]["identity"] >= args.min_identity and match[0]["qcov"] >= args.min_coverage:
            chosen, status = match[0], "anchor_listed_id"
        else:
            status = "listed_id_not_matched"
    if chosen is None:
        full = [h for h in candidates if h["identity"] >= args.min_identity and h["qcov"] >= args.min_coverage
                and h["scov"] >= args.min_coverage]
        in_own = [h for h in full if portal_of(h["subject"]) in own]
        if in_own:
            chosen = max(in_own, key=lambda h: (h["identity"], h["bits"]))
            status = "anchor" if status == "no_hit" else status + ";anchor"
        elif full:
            chosen = max(full, key=lambda h: (h["identity"], h["bits"]))
            status = ("near_identical_other_genome" if own else "near_identical_species_not_in_cohort") + \
                ("" if status == "no_hit" else ";" + status)
        elif candidates:
            chosen = candidates[0]
            status = ("best_hit_only" if own else "best_hit_only_species_not_in_cohort") + \
                ("" if status == "no_hit" else ";" + status)
    is_anchor = status.startswith("anchor") or ";anchor" in status
    if chosen and is_anchor:
        anchor_of.setdefault(chosen["subject"], []).append(sid)
    anchor_rows.append([sid, ref["source"], ref["species"], ref["protein_name"], ref["role"], ref["substrates"],
                        ",".join(sorted(own)), status, chosen["subject"] if chosen else "",
                        round(chosen["identity"], 1) if chosen else "", round(chosen["qcov"], 3) if chosen else "",
                        round(chosen["scov"], 3) if chosen else "",
                        orthogroup.get(chosen["subject"], "") if chosen else ""])

# 2. PF00083 domains: HMM coverage from the union of the domains' HMM coordinates


def read_domtbl(path):
    proteins = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.split()
            p = proteins.setdefault(x[0], dict(length=int(x[2]), hmm_length=int(x[5]), evalue=float(x[6]),
                                               score=float(x[7]), spans=[]))
            p["spans"].append((int(x[15]), int(x[16])))
    for p in proteins.values():
        covered = set()
        for a, b in p["spans"]:
            covered.update(range(a, b + 1))
        p["coverage"] = len(covered) / p["hmm_length"]
        p["domains"] = len(p["spans"])
    return proteins


cohort_hits = read_domtbl(inputs["pf00083_cohort.domtbl"])
reference_hits = read_domtbl(inputs["pf00083_references.domtbl"])
for row in anchor_rows:
    ref_hit = reference_hits.get(row[0])
    row.append(round(ref_hit["coverage"], 3) if ref_hit else "")
full_length = {p for p, h in cohort_hits.items() if h["coverage"] >= args.min_hmm_coverage}

# 3. Orthogroups holding a sugar porter or an anchor
ref_by_id = {r["structure_id"]: r for r in references}
og_rows = []
for og in sorted({orthogroup[p] for p in cohort_hits} | {orthogroup[a] for a in anchor_of if a in orthogroup}):
    genes = members[og]
    species = {portal_of(g) for g in genes}
    phyla = sorted({taxonomy[s]["phylum"] for s in species if s in taxonomy})
    n_full = sum(g in full_length for g in genes)
    n_partial = sum(g in cohort_hits and g not in full_length for g in genes)
    anchored = [sid for g in genes for sid in anchor_of.get(g, [])]
    substrates = sorted({ref_by_id[sid]["substrates"] for sid in anchored if ref_by_id[sid]["source"] == "characterized"
                         and ref_by_id[sid]["role"] == "transporter" and ref_by_id[sid]["substrates"]})
    og_rows.append([og, len(genes), len(species), n_full, n_partial, round(n_full / len(genes), 3), len(anchored),
                    ",".join(anchored), " | ".join(substrates), ",".join(phyla)])
og_rows.sort(key=lambda r: (-r[3], r[0]))

os.makedirs(OUT_DIR, exist_ok=True)


def write(name, header, rows):
    with open(os.path.join(OUT_DIR, name), "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


write("anchors.tsv", ["structure_id", "source", "species", "protein_name", "role", "substrates", "own_portals",
                      "status", "cohort_protein", "identity", "query_coverage", "subject_coverage", "orthogroup",
                      "reference_pf00083_coverage"], anchor_rows)
write("pf00083_hits.tsv", ["clean_id", "portal", "length", "hmm_coverage", "domains", "score", "evalue", "partial",
                           "orthogroup", "anchor_of"],
      [[p, portal_of(p), h["length"], round(h["coverage"], 3), h["domains"], h["score"], h["evalue"],
        p not in full_length, orthogroup.get(p, ""), ",".join(anchor_of.get(p, []))]
       for p, h in sorted(cohort_hits.items())])
per_genome = []
for portal in sorted(taxonomy):
    n_full = sum(1 for p in full_length if portal_of(p) == portal)
    n_partial = sum(1 for p in cohort_hits if portal_of(p) == portal and p not in full_length)
    per_genome.append([portal, taxonomy[portal]["species"], taxonomy[portal]["phylum"], taxonomy[portal]["order"],
                       n_full, n_partial])
write("pf00083_per_genome.tsv", ["portal", "species", "phylum", "order", "sugar_porters", "partial"], per_genome)
write("orthogroup_summary.tsv", ["orthogroup", "genes", "species", "sugar_porters", "partial_sugar_porters",
                                 "share_sugar_porters", "anchors", "anchored_references", "anchor_substrates",
                                 "phyla"], og_rows)
write("st_proteins.tsv", ["clean_id", "portal", "orthogroup", "hmm_coverage", "anchor_of"],
      [[p, portal_of(p), orthogroup.get(p, ""), round(cohort_hits[p]["coverage"], 3), ",".join(anchor_of.get(p, []))]
       for p in sorted(full_length)])

# Report
statuses = {}
for row in anchor_rows:
    statuses[row[7]] = statuses.get(row[7], 0) + 1
print("Anchors: " + ", ".join(f"{k} {v}" for k, v in sorted(statuses.items(), key=lambda x: -x[1])))
shared = {a: s for a, s in anchor_of.items() if len(s) > 1}
if shared:
    print(f"  {len(shared)} cohort proteins anchor more than one reference, e.g. "
          + "; ".join(f"{a}: {','.join(s)}" for a, s in list(shared.items())[:3]))
no_pfam = [r[0] for r in anchor_rows if r[13] == ""]
print(f"References without PF00083: {len(no_pfam)} ({', '.join(no_pfam[:8])}{', ...' if len(no_pfam) > 8 else ''})")
counts = sorted(r[4] for r in per_genome)
print(f"Sugar porters: {len(full_length)} full-length (HMM coverage >= {args.min_hmm_coverage:.0%}) and "
      f"{len(cohort_hits) - len(full_length)} partial, in {len(taxonomy)} genomes "
      f"(median {counts[len(counts) // 2]} per genome, range {counts[0]}-{counts[-1]})")
in_ogs = [r for r in og_rows if r[3] > 0]
singletons = sum(1 for p in full_length if len(members[orthogroup[p]]) == 1)
top = sum(r[3] for r in in_ogs[:10])
print(f"Orthogroups with a full-length sugar porter: {len(in_ogs)} ({singletons} sugar porters unassigned); "
      f"the 10 largest hold {top} ({top / len(full_length):.0%}); "
      f"{sum(1 for r in in_ogs if r[6] == 0)} have no anchored reference")
print(f"Results in {OUT_DIR}/")
