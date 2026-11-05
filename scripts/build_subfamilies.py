#!/usr/bin/env python3
# build_subfamilies.py
# Description: Builds the tip list and sequences of each subfamily tree of a gene family, from the
#              reviewed table families/<family>/subfamilies.tsv, and writes a report to review before
#              alignment. One row per set of tips, with these columns:
#                subfamily   name of the subfamily
#                set         'ingroup' (rows add up), or the name of an outgroup lineage
#                kind        node         the clade with the number in 'definition', from the node table
#                                         of make_family_figure.py for the tree named in 'tree'
#                            clade        the most recent common ancestor of the tips in 'definition',
#                                         in the tree named in 'tree'
#                            orthogroups  all tips of the orthogroups in 'definition'
#                            tips         the tips in 'definition'
#                            uniprot      proteins outside the family tree, as '<tip_id>=<accession>'
#                tree        'ml' or 'fasttree' (kinds node and clade)
#                definition  comma-separated, as above
#                sample      'all', or how many tips to take from an outgroup lineage (kinds node and clade)
#                notes       free text
#              Whole-family trees are rooted on the outgroups with root = TRUE in outgroups.tsv, so
#              that a clade means the same thing in both. A lineage is sampled by cutting its clade
#              into as many parts as tips wanted, always splitting the largest part. Each part gives
#              its characterized protein with the shortest terminal branch, else its member with the
#              shortest terminal branch among those with HMM coverage >= 0.8, against long-branch
#              attraction and fragments. The report does not judge the sets; a person does.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-11-05
# Version: 1.0
#
# Usage (from the project folder; standard library only):
#   python3 scripts/build_subfamilies.py <family> [--lists-only]
#   --lists-only  Tip lists and report without sequences, where tree/sequences.faa is not available.
#                 The full run needs internet for kind = uniprot.
#
# Reads, inside <OUT_DIR> (family.conf): tree/tips.tsv, tree/ml.treefile, tree/fasttree.nwk,
# tree/sequences.faa, members.tsv, figures/<family>_<tree>_nodes.tsv (kind = node); and families/<family>/:
# subfamilies.tsv, outgroups.tsv, references.tsv.
#
# Output, inside <OUT_DIR>/subfamilies:
#   report.md             per subfamily: sets, size, orthogroups, anchors, sampled outgroups; tips in
#                         no subfamily
#   <subfamily>/tips.tsv  tip_id, set, role, kind, portal, orthogroup, references, status
#   <subfamily>/sequences.faa  unaligned, full-length (not with --lists-only)

import argparse
import collections
import hashlib
import os
import sys
import time
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_config import fail, portal_of, read_config, read_fasta, read_tsv, write_tsv  # noqa: E402
from newick_tree import read_tree  # noqa: E402

TREE_FILES = {"ml": "ml.treefile", "fasttree": "fasttree.nwk"}
KINDS = ("node", "clade", "orthogroups", "tips", "uniprot")
TREE_KINDS = ("node", "clade")
MIN_SAMPLE_COVERAGE = 0.8
SIZE_RANGE = (100, 2500)  # the range the plan allows a subfamily tree; outside it the report warns

parser = argparse.ArgumentParser(description="Tip lists and sequences of the subfamily trees of a gene family")
parser.add_argument("family")
parser.add_argument("--lists-only", action="store_true")
args = parser.parse_args()
config = read_config(args.family)
out_dir = config["OUT_DIR"]
tree_dir = os.path.join(out_dir, "tree")
sub_dir = os.path.join(out_dir, "subfamilies")

rows = read_tsv(os.path.join(config["FAMILY_DIR"], "subfamilies.tsv"))
family_tips = {t["tip_id"]: t for t in read_tsv(os.path.join(tree_dir, "tips.tsv"))}
coverage = {m["clean_id"]: float(m["hmm_coverage"]) for m in read_tsv(os.path.join(out_dir, "members.tsv"))}
references = {r["reference_id"]: r for r in read_tsv(os.path.join(config["FAMILY_DIR"], "references.tsv"))}
root_tips = [o["outgroup_id"] for o in read_tsv(os.path.join(config["FAMILY_DIR"], "outgroups.tsv"))
             if o["root"] == "TRUE"]
if not root_tips:
    fail("no outgroup has root = TRUE in outgroups.tsv, so the trees cannot be rooted")


def characterized_names(tip):
    names = []
    for ref in family_tips[tip]["references"].split(","):
        if ref in references and references[ref]["evidence"] == "characterized":
            names.append(references[ref]["protein_name"])
    return names


trees = {}


def get_tree(name):
    """Whole-family tree, rooted on the family's root outgroups."""
    if name not in TREE_FILES:
        fail(f"tree must be one of {', '.join(TREE_FILES)}, not '{name}'")
    if name not in trees:
        path = os.path.join(tree_dir, TREE_FILES[name])
        if not os.path.isfile(path):
            fail(f"'{path}' not found")
        tree = read_tree(path)
        missing = [t for t in root_tips if t not in tree.tips]
        if missing:
            fail(f"root outgroups not in {path}: {', '.join(missing)}")
        tree.root_on_tip(root_tips[0])
        rest = [t for t in tree.tips if t not in root_tips]
        if len(tree.tips_under(tree.mrca(rest))) != len(rest):
            fail(f"the root outgroups {', '.join(root_tips)} do not form a clade in {path}")
        trees[name] = tree
    return trees[name]


def sample_clade(tree, node, wanted, excluded):
    """Cuts a clade into parts and picks one tip per part; returns (tip, why) pairs."""
    sizes = {}
    stack = [node]
    order = []
    while stack:
        current = stack.pop()
        order.append(current)
        stack.extend(current.children)
    for current in reversed(order):
        sizes[id(current)] = sum(sizes[id(c)] for c in current.children) if current.children else 1
    parts = [node]
    while len(parts) < wanted:
        splittable = [p for p in parts if p.children]
        if not splittable:
            break
        largest = max(splittable, key=lambda p: (sizes[id(p)], tree.tips_under(p)[0]))
        parts.remove(largest)
        parts.extend(largest.children)
    picks = []
    for part in parts:
        candidates = [t for t in tree.tips_under(part) if t not in excluded]
        if not candidates:
            continue
        length = lambda t: (tree.tips[t].length, t)  # noqa: E731
        characterized = [t for t in candidates if characterized_names(t)]
        if characterized:
            picks.append((min(characterized, key=length), "characterized"))
            continue
        complete = [t for t in candidates if coverage.get(t, 0) >= MIN_SAMPLE_COVERAGE]
        if complete:
            picks.append((min(complete, key=length), "shortest branch"))
        else:  # only fragments: the most complete one
            picks.append((max(candidates, key=lambda t: (coverage.get(t, 0), t)), "most complete"))
    return picks


node_tables = {}


def node_table(tree_name):
    """Node table of make_family_figure.py for a tree; stops if it belongs to another tree or root."""
    if tree_name not in node_tables:
        path = os.path.join(out_dir, "figures", f"{args.family}_{tree_name}_nodes.tsv")
        if not os.path.isfile(path):
            fail(f"'{path}' not found; run scripts/make_family_figure.py {args.family} --tree {tree_name}")
        table = {r["node"]: r for r in read_tsv(path)}
        with open(os.path.join(tree_dir, TREE_FILES[tree_name]), "rb") as f:
            md5 = hashlib.md5(f.read()).hexdigest()
        for r in table.values():
            if r["tree_md5"] != md5:
                fail(f"'{path}' was made from another version of the tree, so its numbers may mean other clades; "
                     f"run scripts/make_family_figure.py {args.family} --tree {tree_name} and choose the nodes again")
            if r["root"] != ",".join(root_tips):
                fail(f"'{path}' was made with the root '{r['root']}', not {','.join(root_tips)}; "
                     f"run scripts/make_family_figure.py {args.family} --tree {tree_name} and choose the nodes again")
            break
        node_tables[tree_name] = table
    return node_tables[tree_name]


def clade_tips(row):
    """The tips whose most recent common ancestor defines a row of kind node or clade."""
    if row["kind"] == "clade":
        return [x.strip() for x in row["definition"].split(",") if x.strip()]
    number = row["definition"].strip()
    table = node_table(row["tree"])
    if number not in table:
        fail(f"{row['subfamily']}/{row['set']}: node '{number}' is not in the node table of the {row['tree']} tree")
    tree = get_tree(row["tree"])
    entry = table[number]
    if len(tree.tips_under(tree.mrca([entry["tip_a"], entry["tip_b"]]))) != int(entry["tips"]):
        fail(f"{row['subfamily']}/{row['set']}: node {number} does not match the tree; rerun make_family_figure.py")
    return [entry["tip_a"], entry["tip_b"]]


def check_no_root(found, row):
    """A clade holding a root outgroup is not a clade of the rooted tree."""
    inside = [t for t in root_tips if t in found]
    if inside:
        fail(f"{row['subfamily']}/{row['set']}: this clade includes the root outgroup {', '.join(inside)}")


def tips_of_row(row):
    """Tips of a row of the table, as a list (kind = uniprot is handled by the caller)."""
    names = clade_tips(row) if row["kind"] in TREE_KINDS else [x.strip() for x in row["definition"].split(",") if x.strip()]
    if not names:
        fail(f"{row['subfamily']}/{row['set']}: empty definition")
    if row["kind"] == "orthogroups":
        found = [t for t, info in family_tips.items() if info["orthogroup"] in names]
        if not found:
            fail(f"{row['subfamily']}/{row['set']}: no tips in {row['definition']}")
        return found
    unknown = [n for n in names if n not in family_tips]
    if unknown:
        fail(f"{row['subfamily']}/{row['set']}: not tips of the family tree: {', '.join(unknown)}")
    if row["kind"] == "tips":
        return names
    tree = get_tree(row["tree"])
    found = tree.tips_under(tree.mrca(names))
    check_no_root(found, row)
    return found


def top(items, n=5):
    return ", ".join(f"{name} {count}" for name, count in collections.Counter(items).most_common(n))


# 1. Sets of each subfamily
subfamilies = collections.OrderedDict()
for row in rows:
    if row["kind"] not in KINDS:
        fail(f"{row['subfamily']}/{row['set']}: kind must be one of {', '.join(KINDS)}")
    if row["sample"] != "all" and (row["kind"] not in TREE_KINDS or not row["sample"].isdigit() or int(row["sample"]) < 1):
        fail(f"{row['subfamily']}/{row['set']}: sample must be 'all', or a positive number for a node or clade")
    subfamilies.setdefault(row["subfamily"], []).append(row)

results = collections.OrderedDict()
for name, sub_rows in subfamilies.items():
    if not any(r["set"] == "ingroup" for r in sub_rows):
        fail(f"{name}: no ingroup row")
    ingroup, row_notes = {}, []
    for row in sub_rows:
        if row["set"] != "ingroup":
            continue
        if row["sample"] != "all":
            fail(f"{name}: an ingroup takes sample = all")
        found = tips_of_row(row)
        ingroup.update(dict.fromkeys(found))
        support = ""
        if row["kind"] in TREE_KINDS:
            support = get_tree(row["tree"]).mrca(clade_tips(row)).label
        row_notes.append((row, len(found), support))
    lineages = collections.OrderedDict()
    external = []
    for row in sub_rows:
        lineage = row["set"]
        if lineage == "ingroup":
            continue
        if row["kind"] == "uniprot":
            for item in row["definition"].split(","):
                tip_id, _, accession = item.strip().partition("=")
                if not accession:
                    fail(f"{name}/{lineage}: uniprot items must look like '<tip_id>=<accession>'")
                external.append((tip_id, accession, lineage))
            lineages.setdefault(lineage, [])
            continue
        if row["kind"] in TREE_KINDS and row["sample"] != "all":
            tree = get_tree(row["tree"])
            node = tree.mrca(clade_tips(row))
            check_no_root(tree.tips_under(node), row)
            picks = sample_clade(tree, node, int(row["sample"]), set(ingroup))
            row_notes.append((row, len(tree.tips_under(node)), node.label))
        else:
            found = [t for t in tips_of_row(row) if t not in ingroup]
            picks = [(t, "listed") for t in found]
            row_notes.append((row, len(found), ""))
        lineages.setdefault(lineage, []).extend(picks)
    results[name] = {"rows": row_notes, "ingroup": list(ingroup), "lineages": lineages, "external": external}

# 2. Sequences
all_sequences = {}
if not args.lists_only:
    path = os.path.join(tree_dir, "sequences.faa")
    if not os.path.isfile(path):
        fail(f"'{path}' not found; run this where the family's sequences are, or use --lists-only")
    all_sequences = read_fasta(path)


def fetch_uniprot(accession):
    for attempt in range(3):
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
            with urllib.request.urlopen(url, timeout=60) as response:
                text = response.read().decode()
            break
        except OSError:
            time.sleep(5 * (attempt + 1))
    else:
        fail(f"could not download {accession}")
    header, *lines = text.strip().splitlines()
    if f"|{accession}|" not in header:
        fail(f"UniProt returned {header[:60]} for {accession}")
    return "".join(lines)


by_sequence = {}
if all_sequences:
    for tip, sequence in all_sequences.items():
        by_sequence.setdefault(sequence, tip)

# 3. Tip lists and sequences
os.makedirs(sub_dir, exist_ok=True)
for name, result in results.items():
    folder = os.path.join(sub_dir, name)
    os.makedirs(folder, exist_ok=True)
    records = [(t, "ingroup", "ingroup") for t in result["ingroup"]]
    for lineage, picks in result["lineages"].items():
        records.extend((t, lineage, "outgroup") for t, _ in picks)
    extra_sequences = {}
    result["uniprot_status"] = []
    for tip_id, accession, lineage in result["external"]:
        if args.lists_only:
            result["uniprot_status"].append((tip_id, accession, lineage, "not checked (--lists-only)"))
            records.append((tip_id, lineage, "outgroup"))
            continue
        sequence = fetch_uniprot(accession)
        if sequence in by_sequence:  # already in the family: keep the member, under its own ID
            result["uniprot_status"].append((tip_id, accession, lineage, f"is family tip {by_sequence[sequence]}"))
            records.append((by_sequence[sequence], lineage, "outgroup"))
        else:
            result["uniprot_status"].append((tip_id, accession, lineage, "not in the family"))
            extra_sequences[tip_id] = sequence
            records.append((tip_id, lineage, "outgroup"))
    seen, table = set(), []
    for tip, set_name, role in records:
        if tip in seen:
            result.setdefault("repeated", []).append(tip)
            continue
        seen.add(tip)
        info = family_tips.get(tip)
        if info:
            table.append([tip, set_name, role, info["kind"], info["portal"], info["orthogroup"], info["references"],
                          info["status"]])
        else:
            table.append([tip, set_name, role, "uniprot", portal_of(tip), "", "", ""])
    write_tsv(os.path.join(folder, "tips.tsv"),
              ["tip_id", "set", "role", "kind", "portal", "orthogroup", "references", "status"], table)
    if not args.lists_only:
        with open(os.path.join(folder, "sequences.faa"), "w") as f:
            for tip, *_ in table:
                sequence = extra_sequences.get(tip) or all_sequences.get(tip)
                if not sequence:
                    fail(f"{name}: no sequence for {tip}")
                f.write(f">{tip}\n{sequence}\n")
    result["table"] = table

# 4. Report
covered = set()
for result in results.values():
    covered.update(result["ingroup"])
lines = [f"# Subfamilies of {config['NAME']}", "",
         "Written by `scripts/build_subfamilies.py`" + (" with `--lists-only`." if args.lists_only else "."), "",
         "| Subfamily | Ingroup tips | Characterized | *D. squalens* | Outgroup tips | Main orthogroups |",
         "|---|---|---|---|---|---|"]
for name, result in results.items():
    ingroup = result["ingroup"]
    lines.append(f"| {name} | {len(ingroup)} | {sum(bool(characterized_names(t)) for t in ingroup)} "
                 f"| {sum(t.startswith('Dicsqu') for t in ingroup)} | {len(result['table']) - len(ingroup)} "
                 f"| {top([family_tips[t]['orthogroup'] or '(none)' for t in ingroup], 3)} |")
lines.append("")
for name, result in results.items():
    ingroup = result["ingroup"]
    lines += [f"## {name}", ""]
    if not SIZE_RANGE[0] <= len(ingroup) <= SIZE_RANGE[1]:
        lines += [f"**Warning:** {len(ingroup)} ingroup tips, outside the planned {SIZE_RANGE[0]}–{SIZE_RANGE[1]}.", ""]
    others = [(o, len(set(ingroup) & set(r["ingroup"]))) for o, r in results.items() if o != name]
    shared = [f"{o} ({n})" for o, n in others if n]
    if shared:
        lines += [f"**Ingroup tips also in:** {', '.join(shared)}.", ""]
    lines += ["| Set | Kind | Tree | Tips in the set | Support | Sample |", "|---|---|---|---|---|---|"]
    for row, size, support in result["rows"]:
        kind = f"node {row['definition']}" if row["kind"] == "node" else row["kind"]
        lines.append(f"| {row['set']} | {kind} | {row['tree']} | {size} | {support} | {row['sample']} |")
    names = sorted({n for t in ingroup for n in characterized_names(t)})
    lines += ["", f"Orthogroups: {top([family_tips[t]['orthogroup'] or '(none)' for t in ingroup], 6)}.", "",
              f"Characterized ({len(names)}): {', '.join(names) if names else 'none'}.", ""]
    for lineage, picks in result["lineages"].items():
        if picks:
            lines.append(f"Outgroup `{lineage}`: " + ", ".join(
                f"{t}" + (f" ({'/'.join(characterized_names(t))})" if characterized_names(t) else "")
                for t, _ in picks) + ".")
            lines.append("")
    for tip_id, accession, lineage, status in result["uniprot_status"]:
        lines.append(f"- `{lineage}`: {tip_id} (UniProt {accession}): {status}")
    if result["uniprot_status"]:
        lines.append("")
    if result.get("repeated"):
        lines += [f"**Listed twice, kept once:** {', '.join(result['repeated'])}.", ""]

left = [t for t in family_tips if t not in covered]
left_characterized = sorted(n for t in left for n in characterized_names(t))
lines += ["## Tips in no subfamily", "",
          f"{len(left)} of {len(family_tips)} tips are in no ingroup, "
          f"{sum(t.startswith('Dicsqu') for t in left)} of them *D. squalens*.", "",
          f"Orthogroups: {top([family_tips[t]['orthogroup'] or '(none)' for t in left], 12)}.", "",
          f"Characterized: {', '.join(left_characterized) if left_characterized else 'none'}.", ""]
report = os.path.join(sub_dir, "report.md")
with open(report, "w") as f:
    f.write("\n".join(lines))
print(f"{len(results)} subfamilies; {len(covered)} of {len(family_tips)} tips in an ingroup. Review {report}")
if args.lists_only:
    print("Tip lists only. Run without --lists-only where tree/sequences.faa is, before aligning.")
