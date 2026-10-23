#!/usr/bin/env python3
# summarize_family.py
# Description: Turns the searches of search_family.sbatch into a gene family's anchors, orphans,
#              members and orthogroups.
#              1. Anchors: the cohort copy of each reference protein (references.tsv, use = TRUE),
#                 from its DIAMOND hits, by the first rule that applies:
#                   anchor_listed_id           the clean ID listed for it (cohort_id), with >= 98%
#                                              identity over >= 90% of the reference;
#                   anchor                     in its own genome, >= 98% identity over >= 90% of both
#                                              proteins;
#                   anchor_gene_model_differs  the same, but over >= 90% of the shorter protein only;
#                   probable_anchor            in its own genome, 95-98% identity over >= 90% of the
#                                              shorter protein, and the cohort protein's closest
#                                              reference is this one;
#                   orphan                     otherwise: species not in the cohort, or no copy in its
#                                              genome. Orphans enter the trees with their own sequence
#                                              as '<short_name>-<protein name>', where short_name is a
#                                              species-only pseudo-portal (families/<family>/orphans.tsv
#                                              gives its taxonomy).
#                 The own genome is the reference's cohort_portal, else the cohort genome(s) of its
#                 species (also under the names given as '(synonym ...)' or '(anamorph ...)').
#              2. Members: cohort proteins whose domains of one family HMM cover >= MIN_HMM_COVERAGE of
#                 the HMM, at the gathering threshold; hits below that are kept in hmm_hits.tsv.
#              3. Orthogroups (OrthoFinder Orthogroups.tsv, the root hierarchical orthogroups) holding a
#                 member or an anchor.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-10-23
# Version: 1.0
#
# Usage (from the project folder; standard library only): python3 scripts/summarize_family.py <family>
#
# Output, inside OUT_DIR (family.conf):
#   anchors.tsv              per reference: status, cohort protein (or best hit), identity, coverages,
#                            orthogroup, tree ID (the anchor, or the orphan's ID) and HMM coverage
#   hmm_hits.tsv             per cohort protein with a family HMM hit: coverage, score, member or not
#   members_per_genome.tsv   members and sub-threshold hits per cohort genome
#   orthogroup_summary.tsv   per orthogroup with a member or an anchor
#   members.tsv              the members, with orthogroup and anchor status (tree input)
#   orphans_missing.tsv      orphan pseudo-portals not yet described in families/<family>/orphans.tsv

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_config import OG_DIR, fail, portal_of, read_config, read_tsv, taxonomy, write_tsv  # noqa: E402

ANCHOR_IDENTITY, PROBABLE_IDENTITY, MIN_COVERAGE = 98.0, 95.0, 0.9
# Older names in reference tables for species named differently in the cohort taxonomy
SYNONYMS = {"piriformospora indica": "serendipita indica", "postia placenta": "rhodonia placenta",
            "myceliophthora thermophila": "thermothelomyces thermophilus",
            "pichia stipitis": "scheffersomyces stipitis", "magnaporthe oryzae": "pyricularia oryzae",
            "dekkera bruxellensis": "brettanomyces bruxellensis"}

if len(sys.argv) != 2:
    fail("usage: python3 scripts/summarize_family.py <family>")
config = read_config(sys.argv[1])
out, search = config["OUT_DIR"], os.path.join(config["OUT_DIR"], "search")
references = [r for r in read_tsv(os.path.join(config["FAMILY_DIR"], "references.tsv")) if r["use"] == "TRUE"]
orphan_table = {r["pseudo_portal"]: r for r in read_tsv(os.path.join(config["FAMILY_DIR"], "orphans.tsv"),
                                                          required=False)}
tax = taxonomy()


def binomial(name):
    return " ".join(re.sub(r"[^A-Za-z. -]", " ", name).split()[:2]).lower()


portals_of_species = {}
for portal, row in tax.items():
    portals_of_species.setdefault(binomial(row["species"]), set()).add(portal)


def own_portals(ref):
    """The reference's own cohort genome: its cohort_portal if that is a cohort genome, else by species."""
    if ref["cohort_portal"] in tax:
        return {ref["cohort_portal"]}
    if ref["cohort_portal"]:
        print(f"  note: {ref['reference_id']}: cohort_portal {ref['cohort_portal']} is not a cohort genome; "
              f"using its species instead")
    names = set()
    for part in ref["species"].split(" / "):
        names.add(binomial(part.split("(")[0]))
        names |= {binomial(a) for a in re.findall(r"\((?:synonym|anamorph)\s+([^)]+)\)", part)}
    names |= {SYNONYMS[n] for n in names if n in SYNONYMS}
    return set().union(*(portals_of_species.get(n, set()) for n in names))


def safe(text):
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]+", "_", text)).strip("_")


# Orthogroup of every cohort gene; unassigned genes are single-gene orthogroups
orthogroup, members_of = {}, {}
for name in ("Orthogroups.tsv", "Orthogroups_UnassignedGenes.tsv"):
    path = os.path.join(OG_DIR, name)
    if not os.path.isfile(path):
        fail(f"'{path}' not found; copy the OrthoFinder Orthogroups folder there")
    with open(path) as f:
        f.readline()
        for line in f:
            fields = line.rstrip("\r\n").split("\t")
            genes = [g for cell in fields[1:] for g in cell.split(", ") if g]
            members_of[fields[0]] = genes
            for gene in genes:
                orthogroup[gene] = fields[0]

# 1. Anchors and orphans
diamond = os.path.join(search, "diamond_references.tsv")
if not os.path.isfile(diamond):
    fail(f"'{diamond}' not found; run scripts/search_family.sbatch {config['FAMILY']}")
hits = {}
best_reference = {}  # cohort protein -> (bitscore, reference): its closest reference
with open(diamond) as f:
    for line in f:
        q, s, pident, length, qlen, slen, qs, qe, ss, se, evalue, bits = line.rstrip("\n").split("\t")
        qcov, scov = (int(qe) - int(qs) + 1) / int(qlen), (int(se) - int(ss) + 1) / int(slen)
        short, long_ = (qcov, scov) if int(qlen) <= int(slen) else (scov, qcov)
        hits.setdefault(q, []).append(dict(subject=s, identity=float(pident), qcov=qcov, scov=scov,
                                           short_cov=short, long_cov=long_, bits=float(bits)))
        if float(bits) > best_reference.get(s, (0, ""))[0]:
            best_reference[s] = (float(bits), q)

anchor_rows, anchor_of, status_of, orphans = [], {}, {}, {}
for ref in references:
    rid, own = ref["reference_id"], own_portals(ref)
    candidates = sorted(hits.get(rid, []), key=lambda h: -h["bits"])
    in_own = [h for h in candidates if portal_of(h["subject"]) in own]
    chosen, status = None, None
    listed = [h for h in candidates if h["subject"] == ref["cohort_id"]] if ref["cohort_id"] else []
    if listed and listed[0]["identity"] >= ANCHOR_IDENTITY and listed[0]["qcov"] >= MIN_COVERAGE:
        chosen, status = listed[0], "anchor_listed_id"
    if chosen is None:
        strict = [h for h in in_own if h["identity"] >= ANCHOR_IDENTITY and h["short_cov"] >= MIN_COVERAGE]
        if strict:
            chosen = max(strict, key=lambda h: (h["long_cov"] >= MIN_COVERAGE, h["identity"], h["bits"]))
            status = "anchor" if chosen["long_cov"] >= MIN_COVERAGE else "anchor_gene_model_differs"
    if chosen is None:
        probable = [h for h in in_own if PROBABLE_IDENTITY <= h["identity"] < ANCHOR_IDENTITY
                    and h["short_cov"] >= MIN_COVERAGE and best_reference[h["subject"]][1] == rid]
        if probable:
            chosen, status = max(probable, key=lambda h: (h["identity"], h["bits"])), "probable_anchor"
    if chosen is None:
        chosen = candidates[0] if candidates else None
        status = "orphan_no_copy_in_genome" if own else "orphan_species_not_in_cohort"
        pseudo = safe(ref["short_name"])
        if pseudo in tax:
            fail(f"orphan pseudo-portal '{pseudo}' ({rid}) is a cohort portal; change its short_name")
        if pseudo in orphans and binomial(orphans[pseudo]["species"]) != binomial(ref["species"]):
            fail(f"pseudo-portal '{pseudo}' is used for two species: {orphans[pseudo]['species']} and {ref['species']}")
        orphans.setdefault(pseudo, dict(species=ref["species"], ids=[]))
        tree_id = f"{pseudo}-{safe(ref['protein_name'])}"
        if tree_id in {i for o in orphans.values() for i in o["ids"]}:
            fail(f"two orphans would both be named '{tree_id}'")
        orphans[pseudo]["ids"].append(tree_id)
    else:
        anchor_of.setdefault(chosen["subject"], []).append(rid)
        tree_id = chosen["subject"]
    status_of[rid] = status
    anchor_rows.append([rid, ref["evidence"], ref["species"], ref["protein_name"], ref["role"], ref["function"],
                        ",".join(sorted(own)), status, chosen["subject"] if chosen else "",
                        round(chosen["identity"], 1) if chosen else "", round(chosen["qcov"], 3) if chosen else "",
                        round(chosen["scov"], 3) if chosen else "",
                        orthogroup.get(chosen["subject"], "") if chosen else "", tree_id])

# 2. HMM hits: coverage from the union of each HMM's domain coordinates


def read_domtbl(path):
    proteins = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.split()
            p = proteins.setdefault(x[0], dict(length=int(x[2]), hmm=x[4], hmm_length=int(x[5]),
                                               evalue=float(x[6]), score=float(x[7]), spans=[]))
            p["spans"].append((int(x[15]), int(x[16])))
    for p in proteins.values():
        covered = set()
        for a, b in p["spans"]:
            covered.update(range(a, b + 1))
        p["coverage"] = len(covered) / p["hmm_length"]
    return proteins


cohort_hits, reference_cov = {}, {}
for hmm in config["HMMS"]:
    for target in ("cohort", "references"):
        path = os.path.join(search, f"{hmm.lower()}_{target}.domtbl")
        if not os.path.isfile(path):
            fail(f"'{path}' not found; run scripts/search_family.sbatch {config['FAMILY']}")
        for protein, hit in read_domtbl(path).items():
            store = cohort_hits if target == "cohort" else reference_cov
            # With several HMMs, a protein keeps its best-covering one
            if target == "cohort":
                if protein not in store or hit["coverage"] > store[protein]["coverage"]:
                    store[protein] = hit
            else:
                store[protein] = max(store.get(protein, 0.0), hit["coverage"])
for row in anchor_rows:
    row.append(round(reference_cov[row[0]], 3) if row[0] in reference_cov else "")
members = {p for p, h in cohort_hits.items() if h["coverage"] >= config["MIN_HMM_COVERAGE"]}

# 3. Orthogroups holding a member or an anchor
ref_by_id = {r["reference_id"]: r for r in references}
og_rows = []
for og in sorted({orthogroup[p] for p in cohort_hits if p in orthogroup} | {orthogroup[a] for a in anchor_of}):
    genes = members_of[og]
    species = {portal_of(g) for g in genes}
    n_members = sum(g in members for g in genes)
    n_below = sum(g in cohort_hits and g not in members for g in genes)
    anchored = [rid for g in genes for rid in anchor_of.get(g, [])]
    functions = sorted({ref_by_id[r]["function"] for r in anchored if ref_by_id[r]["evidence"] == "characterized"
                        and ref_by_id[r]["function"]})
    og_rows.append([og, len(genes), len(species), n_members, n_below, round(n_members / len(genes), 3), len(anchored),
                    ",".join(anchored), " | ".join(functions),
                    ",".join(sorted({tax[s]["phylum"] for s in species if s in tax}))])
og_rows.sort(key=lambda r: (-r[3], r[0]))

write_tsv(os.path.join(out, "anchors.tsv"),
          ["reference_id", "evidence", "species", "protein_name", "role", "function", "own_portals", "status",
           "cohort_protein", "identity", "query_coverage", "subject_coverage", "orthogroup", "tree_id",
           "reference_hmm_coverage"], anchor_rows)
write_tsv(os.path.join(out, "hmm_hits.tsv"),
          ["clean_id", "portal", "length", "hmm", "hmm_coverage", "domains", "score", "evalue", "member",
           "orthogroup", "anchor_of"],
          [[p, portal_of(p), h["length"], h["hmm"], round(h["coverage"], 3), len(h["spans"]), h["score"], h["evalue"],
            p in members, orthogroup.get(p, ""), ",".join(anchor_of.get(p, []))] for p, h in sorted(cohort_hits.items())])
per_genome = [[p, tax[p]["species"], tax[p]["phylum"], tax[p]["order"], sum(portal_of(m) == p for m in members),
               sum(portal_of(h) == p for h in cohort_hits if h not in members)] for p in sorted(tax)]
write_tsv(os.path.join(out, "members_per_genome.tsv"),
          ["portal", "species", "phylum", "order", "members", "below_coverage"], per_genome)
write_tsv(os.path.join(out, "orthogroup_summary.tsv"),
          ["orthogroup", "genes", "species", "members", "below_coverage", "share_members", "anchors",
           "anchored_references", "anchor_functions", "phyla"], og_rows)
write_tsv(os.path.join(out, "members.tsv"),
          ["clean_id", "portal", "orthogroup", "hmm", "hmm_coverage", "anchor_of", "anchor_status"],
          [[p, portal_of(p), orthogroup.get(p, ""), cohort_hits[p]["hmm"], round(cohort_hits[p]["coverage"], 3),
            ",".join(anchor_of.get(p, [])), ",".join(sorted({status_of[r] for r in anchor_of.get(p, [])}))]
           for p in sorted(members)])
missing = [[p, o["species"], ",".join(o["ids"])] for p, o in sorted(orphans.items()) if p not in orphan_table]
write_tsv(os.path.join(out, "orphans_missing.tsv"), ["pseudo_portal", "species", "tree_ids"], missing)

# Report
counts = {}
for s in status_of.values():
    counts[s] = counts.get(s, 0) + 1
print("References: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items(), key=lambda x: -x[1])))
shared = {a: r for a, r in anchor_of.items() if len(r) > 1}
if shared:
    print(f"  {len(shared)} cohort proteins anchor several references: "
          + "; ".join(f"{a}: {','.join(r)}" for a, r in shared.items()))
not_member = [r[0] for r in anchor_rows if float(r[14] or 0) < config["MIN_HMM_COVERAGE"]]
print(f"References below the membership rule themselves: {', '.join(not_member) or 'none'}")
n = sorted(r[4] for r in per_genome)
print(f"Members: {len(members)} ({len(cohort_hits) - len(members)} hits below {config['MIN_HMM_COVERAGE']:.0%} "
      f"coverage), median {n[len(n) // 2]} per genome (range {n[0]}-{n[-1]}), in "
      f"{sum(1 for r in og_rows if r[3] > 0)} orthogroups")
print(f"Orphans: {sum(len(o['ids']) for o in orphans.values())} proteins under {len(orphans)} pseudo-portals"
      + (f"; {len(missing)} pseudo-portals still to describe in families/{config['FAMILY']}/orphans.tsv "
         f"(listed in {out}/orphans_missing.tsv)" if missing else ""))
print(f"Results in {out}/")
