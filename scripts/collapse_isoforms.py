#!/usr/bin/env python3
"""collapse_isoforms.py

Description: Collapses isoforms and redundant gene models of one genome (step A of the proteome
             cleanup). MycoCosm GFFs store every isoform as a separate gene, so models that share
             most of their coding sequence (CDS) on the same strand are grouped into one locus, and
             only the longest protein of each locus is kept. The clean proteome is the input
             proteome filtered by proteinId (nothing is translated).
Author: Victor M. Gonzalez Ramos
Date: 2026-10-15
Version: 1.0

Usage: collapse_isoforms.py --portal P --gff P_genes.gff3.gz --proteome P_proteins.faa \
           --work-dir work/P --out-fasta collapsed/P_proteins.faa [--min-shared-pct 50]
Normally run through collapse_isoforms.sh, which loops over the portal list.

Rules:
- Two models are grouped when, on the same strand, they share at least --min-shared-pct of the
  shorter model's CDS, also through a chain (A with B, B with C). Overlapping pairs below the
  threshold stay separate genes and are listed in unmerged.tsv.
- The longest protein of each locus is kept. Ties are broken by the gene name (alphabetically, so
  the Broad T0 transcript is kept), then by the lowest proteinId.
- Proteins that are in the proteome but not in the GFF are kept as they are.

Writes to --work-dir: loci.tsv, pairs.tsv, removed.tsv, unmerged.tsv and summary_row.tsv, and the
collapsed proteome to --out-fasta. The proteome is only written once all checks pass.
"""

import argparse
import gzip
import os
import sys

import networkx as nx
import numpy as np
import pandas as pd


def fail(message):
    sys.exit(f"ERROR: {message}")


def parse_attributes(column9):
    """'ID=a;Parent=b' -> {'ID': 'a', 'Parent': 'b'}"""
    return dict(field.split("=", 1) for field in column9.strip().split(";") if "=" in field)


def read_gff(path):
    """Reads the gene, mRNA and CDS lines of a MycoCosm GFF3.

    Returns a table of CDS pieces (scaffold, strand, start, end, proteinId) and a dict of
    proteinId -> gene name. Only lines starting with '#' are comments: gene names can contain
    '#' (e.g. fgenesh1_kg.1_#_3_#_transcript), so '#' elsewhere in a line is not a comment.
    """
    names = {}
    mrna_protein = {}
    cds_rows = []
    with gzip.open(path, "rt") as gff:
        for line in gff:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                fail(f"{path}: line with {len(fields)} columns instead of 9: {line[:100]}")
            scaffold, _, feature, start, end, _, strand, _, column9 = fields
            if feature == "gene":
                attributes = parse_attributes(column9)
                if "proteinId" in attributes and "feature_name" in attributes:
                    names[attributes["proteinId"]] = attributes["feature_name"]
            elif feature == "mRNA":
                attributes = parse_attributes(column9)
                if "proteinId" not in attributes:
                    fail(f"{path}: mRNA without proteinId: {attributes.get('ID')}")
                mrna_protein[attributes["ID"]] = attributes["proteinId"]
            elif feature == "CDS":
                # GFF is 1-based and closed; start - 1 makes it 0-based and half-open like BED,
                # so the length of a piece is end - start
                cds_rows.append((scaffold, strand, int(start) - 1, int(end), parse_attributes(column9).get("Parent")))

    if not cds_rows:
        fail(f"{path}: no CDS lines")
    cds = pd.DataFrame(cds_rows, columns=["scaffold", "strand", "start", "end", "mrna"])
    unknown = set(cds["mrna"]) - set(mrna_protein)
    if unknown:
        fail(f"{path}: {len(unknown)} CDS Parent IDs are not mRNAs, e.g. {sorted(unknown)[:3]}")
    cds["protein_id"] = cds["mrna"].map(mrna_protein)
    return cds.drop(columns="mrna"), names


def read_proteome(path):
    """Returns {proteinId: (header line, sequence lines as written, protein length)}, in file order.

    The proteinId is the 3rd '|' field of the header: >jgi|Neucr2|4276|NCU00006T0 ...
    """
    records = {}
    protein_id = None
    with open(path) as fasta:
        for line in fasta:
            if line.startswith(">"):
                fields = line[1:].split("|")
                if len(fields) < 3:
                    fail(f"{path}: header without a proteinId field: {line.strip()}")
                protein_id = fields[2]
                if protein_id in records:
                    fail(f"{path}: proteinId {protein_id} appears twice")
                records[protein_id] = [line, [], 0]
            elif protein_id is not None:
                records[protein_id][1].append(line)
                records[protein_id][2] += len(line.strip())
    return records


def overlapping_pairs(cds):
    """Shared CDS bp for every pair of models with overlapping CDS on the same strand.

    Within each scaffold and strand, CDS pieces are sorted by start, so the pieces that overlap
    piece i are the ones right after it that start before piece i ends. Shared bp are added up
    per pair of models, since two models can overlap with several of their CDS pieces.
    """
    shared = {}
    for _, group in cds.sort_values(["scaffold", "strand", "start"]).groupby(["scaffold", "strand"], sort=False):
        start = group["start"].to_list()
        end = group["end"].to_list()
        ids = group["protein_id"].to_list()
        for i in range(len(start)):
            j = i + 1
            while j < len(start) and start[j] < end[i]:
                if ids[i] != ids[j]:
                    pair = tuple(sorted((ids[i], ids[j]), key=int))
                    shared[pair] = shared.get(pair, 0) + min(end[i], end[j]) - start[j]
                j += 1

    pairs = pd.DataFrame([(a, b, bp) for (a, b), bp in shared.items()],
                         columns=["protein_a", "protein_b", "shared_bp"])
    pairs = pairs.sort_values(["protein_a", "protein_b"], key=lambda column: column.astype(int))

    cds_length = (cds["end"] - cds["start"]).groupby(cds["protein_id"]).sum()
    shorter = np.minimum(pairs["protein_a"].map(cds_length), pairs["protein_b"].map(cds_length))
    pairs["pct_of_shorter_cds"] = 100 * pairs["shared_bp"] / shorter
    return pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--portal", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--proteome", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--min-shared-pct", type=float, default=50)
    args = parser.parse_args()
    os.makedirs(args.work_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.out_fasta) or ".", exist_ok=True)

    cds, names = read_gff(args.gff)
    proteome = read_proteome(args.proteome)

    models = sorted(set(cds["protein_id"]), key=int)
    missing = [p for p in models if p not in proteome]
    if missing:
        fail(f"{len(missing)} proteinIds are in the GFF but not in the proteome, e.g. {missing[:3]}")

    # Pairs of models with overlapping CDS, and loci = connected groups of pairs above the threshold
    pairs = overlapping_pairs(cds)
    merged = pairs[pairs["pct_of_shorter_cds"] >= args.min_shared_pct]
    unmerged = pairs[pairs["pct_of_shorter_cds"] < args.min_shared_pct]
    graph = nx.Graph()
    graph.add_nodes_from(models)
    graph.add_edges_from(zip(merged["protein_a"], merged["protein_b"]))

    # One row per model: its locus (lowest proteinId of the locus) and the model kept for it.
    # Sorting by length (longest first), gene name, then proteinId puts the kept model first.
    rows = []
    for locus in nx.connected_components(graph):
        ranked = sorted(locus, key=lambda p: (-proteome[p][2], names.get(p, ""), int(p)))
        locus_id = min(locus, key=int)
        rows += [(p, locus_id, proteome[p][2], ranked[0]) for p in locus]
    loci = pd.DataFrame(rows, columns=["protein_id", "locus", "length", "kept_protein_id"])
    loci = loci.sort_values(["locus", "protein_id"], key=lambda column: column.astype(int))

    kept = set(loci.loc[loci["protein_id"] == loci["kept_protein_id"], "protein_id"])
    not_in_gff = [p for p in proteome if p not in graph]
    keep = kept | set(not_in_gff)
    removed = loci[loci["protein_id"] != loci["kept_protein_id"]]

    # Intermediate tables, with the portal in front where they are combined across genomes
    pairs.to_csv(f"{args.work_dir}/pairs.tsv", sep="\t", index=False, float_format="%.1f")
    loci.to_csv(f"{args.work_dir}/loci.tsv", sep="\t", index=False)
    pd.DataFrame({
        "portal": args.portal,
        "removed_proteinId": removed["protein_id"],
        "kept_proteinId": removed["kept_protein_id"],
        "removed_length": removed["length"],
        "kept_length": removed["kept_protein_id"].map(lambda p: proteome[p][2]),
    }).to_csv(f"{args.work_dir}/removed.tsv", sep="\t", index=False, header=False)
    unmerged.assign(portal=args.portal)[["portal", "protein_a", "protein_b", "shared_bp", "pct_of_shorter_cds"]].to_csv(
        f"{args.work_dir}/unmerged.tsv", sep="\t", index=False, header=False, float_format="%.1f")

    # Collapsed proteome, in the input order, written to .part until the checks pass
    part = f"{args.out_fasta}.part"
    with open(part, "w") as out:
        for protein_id, (header, sequence, _) in proteome.items():
            if protein_id in keep:
                out.write(header)
                out.writelines(sequence)
    n_out = sum(1 for line in open(part) if line.startswith(">"))
    if n_out != len(keep) or len(proteome) - n_out != len(removed):
        os.remove(part)
        fail(f"counts do not add up (in {len(proteome)}, out {n_out}, kept {len(keep)}, removed {len(removed)})")

    # summary row: portal, proteins in, GFF models, loci, removed, % removed, not in GFF, proteins out,
    # overlapping pairs left unmerged, and the largest % of CDS shared by an unmerged pair
    max_unmerged = unmerged["pct_of_shorter_cds"].max() if len(unmerged) else 0
    summary = [args.portal, len(proteome), len(models), len(kept), len(removed),
               f"{100 * len(removed) / len(models):.2f}", len(not_in_gff), n_out, len(unmerged), f"{max_unmerged:.1f}"]
    with open(f"{args.work_dir}/summary_row.tsv", "w") as out:
        out.write("\t".join(map(str, summary)) + "\n")
    os.replace(part, args.out_fasta)

    print(f"  {args.portal}: {len(models)} models -> {len(kept)} loci ({len(removed)} removed, "
          f"{len(unmerged)} overlapping pairs below {args.min_shared_pct:g}% left unmerged, "
          f"{len(not_in_gff)} proteins not in the GFF)")


if __name__ == "__main__":
    main()
