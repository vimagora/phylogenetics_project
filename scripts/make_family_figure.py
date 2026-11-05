#!/usr/bin/env python3
# make_family_figure.py
# Description: Draws a gene family's tree as a circular figure (SVG), for slides and for review. From
#              the centre outward: the tree; one ring coloured by orthogroup (those with at least 1% of
#              the tips, in the colours of the iTOL files); per family HMM, a ring of the full-sequence
#              bit score and a ring of the HMM's coverage; ticks for the proteins of the portals in
#              MARK_PORTALS (family.conf, optional); and the names of characterized and genetically
#              supported references (black) and of outgroups (red). Names are moved apart along the
#              circle as little as possible, with leader lines to their tips. Nothing is inferred for
#              other tips.
#              The tree is rooted on the outgroups with root = TRUE in outgroups.tsv, or at its midpoint
#              if there are none or they do not form a clade; the caption says which.
#              For choosing subfamilies (step F10), it also writes a node figure and a node table. The
#              node figure has the same rings and names, but its tree shows only the clades of at least
#              --min-node-tips tips, numbered by size, evenly spaced to show how they nest (not branch
#              length), each in a circle coloured by support. Support is a flag, not a filter. The table
#              describes each numbered clade; subfamilies.tsv can name them by number (kind = node).
#              With --subfamily it draws that subfamily's tree instead (step F14), rooted on its outgroups,
#              whose names are red and whose tips without a name show their tip ID, with no node figure.
#              With --subfamilies it also writes a copy of the figure showing where the subfamilies
#              landed: each ingroup is a semitransparent slice behind the tree, starting at its common
#              ancestor and stopping at the rings (so their colours stay true), continued by an arc outside
#              them and named on its largest piece. A group that is not a clade in this tree shows as
#              several pieces, and overlaps as darker areas.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-11-05
# Version: 1.1
#
# Usage (from the project folder; standard library only):
#   python3 scripts/make_family_figure.py <family> [--subfamily NAME | --subfamilies [A,B,...]] [--tree fasttree|ml]
#                                         [--out FILE] [--min-node-tips N]
#   --subfamily      The subfamily to draw, a folder in <OUT_DIR>/subfamilies (default: the whole family).
#   --subfamilies    Also write <figure>_subfamilies.svg, with slices for the listed subfamilies (default: all
#                    in subfamilies.tsv, which need their tips.tsv from build_subfamilies.py).
#   --tree           The tree to draw (default: fasttree).
#   --out            Where to write the figure (default: <OUT_DIR>/figures/<family>[_<subfamily>]_<tree>.svg); the node
#                    figure goes next to it, as <name>_nodes.svg.
#   --min-node-tips  Smallest clade to number (default: 100, the smallest planned subfamily).
#
# Reads, inside <OUT_DIR> (family.conf): tree/tips.tsv, tree/fasttree.nwk or tree/ml.treefile (with
# --subfamily, the same files in subfamilies/<subfamily>), and
# search/<hmm>_cohort.domtbl and search/<hmm>_references.domtbl for each HMM; and from
# families/<family>/: references.tsv, outgroups.tsv (optional).
#
# Output: the figure, the node figure, and <OUT_DIR>/figures/<family>_<tree>_nodes.tsv: node, tips, support,
# main_orthogroups, main_share, named (names in the clade), marked (MARK_PORTALS proteins), tip_a and tip_b
# (two tips whose most recent common ancestor is the clade), tree_md5 and root (what the numbers belong to).

import argparse
import hashlib
import math
import os
import re
import sys
from xml.sax.saxutils import escape

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_config import (fail, orthogroup_colours, read_config, read_domtbl, read_tsv, taxonomy,  # noqa: E402
                           write_tsv)
from newick_tree import read_tree  # noqa: E402

TREES = {"fasttree": ("fasttree.nwk", "FastTree"), "ml": ("ml.treefile", "IQ-TREE")}
LABEL_EVIDENCE = ("characterized", "genetic")
WIDTH = 900.0
OUTER_RADIUS = 435.0            # the circle, names included
BAND = 36.0                     # total width of the rings
RING_GAP = 1.2
TICK_TRACK = 6.0
LEADER = 26.0                   # room between the last ring and the names
FONT_SIZE = 3.4
MIN_FONT_SIZE = 2.4             # below this, the names are too many to read on a slide: a warning is printed
LINE_HEIGHT = 1.3               # space per name along the circle, in font sizes
TIP_QUANTILE = 0.999            # the rings sit just outside this share of the tips
OTHER_ORTHOGROUP = "#d9d9d9"
NO_SCORE = "#ececec"
HEAT_BINS = 32
SCORE_RAMP = ("#ffffcc", "#fd8d3c", "#800026")
COVERAGE_RAMP = ("#f7fcf0", "#7bccc4", "#084081")
ROOT_COLOUR = "#b2182b"
SCORE_QUANTILE = 0.99           # top of each score ring, so that outliers do not wash it out
SLICE_COLOURS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
SLICE_EXTRA_COLOUR = "#8c8c8c"   # subfamilies beyond the colours; each slice carries its name
SLICE_OPACITY = 0.2
SLICE_ARC = 2.0
SLICE_FONT_SIZE = 8.0
NODE_RADIUS = 3.4
NODE_FONT_SIZE = 2.9
SUPPORT_COLOURS = ("#f46d43", "#fee08b", "#66bd63")  # low, middle, high
SUPPORT_CLASSES = {"fasttree": (("< 0.7", "0.7–0.9", "≥ 0.9"), 0.7, 0.9, lambda label: float(label)),
                   "ml": (("UFBoot < 70", "UFBoot 70–95", "UFBoot ≥ 95"), 70, 95,
                          lambda label: float(label.split("/")[-1]))}


def hex_rgb(colour):
    return tuple(int(colour[i:i + 2], 16) for i in (1, 3, 5))


def ramp_colour(ramp, x):
    """Colour at x (0-1) along a ramp of three colours."""
    x = min(max(x, 0.0), 1.0)
    low, high, t = (ramp[0], ramp[1], x * 2) if x < 0.5 else (ramp[1], ramp[2], x * 2 - 1)
    a, b = hex_rgb(low), hex_rgb(high)
    return "#" + "".join(f"{round(p + (q - p) * t):02x}" for p, q in zip(a, b))


def species_prefix(species):
    words = re.sub(r"\(.*?\)", "", species).split()
    return (words[0][0] + words[1][0]).capitalize() if len(words) > 1 else words[0][:2]


def spread(positions, gap, low, high):
    """Positions as close as possible to the given sorted ones, at least gap apart, inside low-high:
    least squares under an order constraint (pool-adjacent-violators), then clamped."""
    blocks = []  # [sum of shifted positions, count]
    for i, p in enumerate(positions):
        blocks.append([p - i * gap, 1])
        while len(blocks) > 1 and blocks[-2][0] / blocks[-2][1] > blocks[-1][0] / blocks[-1][1]:
            total, count = blocks.pop()
            blocks[-1][0] += total
            blocks[-1][1] += count
    shifted = [total / count for total, count in blocks for _ in range(count)]
    n = len(positions)
    return [min(max(s + i * gap, low + i * gap), high - (n - 1 - i) * gap) for i, s in enumerate(shifted)]


parser = argparse.ArgumentParser(description="Circular figure of a gene family's tree")
parser.add_argument("family")
parser.add_argument("--subfamily")
parser.add_argument("--subfamilies", nargs="?", const="all")
parser.add_argument("--tree", choices=TREES, default="fasttree")
parser.add_argument("--out")
parser.add_argument("--min-node-tips", type=int, default=100)
args = parser.parse_args()
if args.subfamily and args.subfamilies:
    fail("--subfamilies adds slices to the whole-family figure, so it cannot go with --subfamily")
config = read_config(args.family)
planned, chosen = [], []
if args.subfamilies:
    planned = list(dict.fromkeys(r["subfamily"] for r in read_tsv(os.path.join(config["FAMILY_DIR"], "subfamilies.tsv"))))
    chosen = planned if args.subfamilies == "all" else args.subfamilies.split(",")
    for name in chosen:
        if name not in planned:
            fail(f"subfamily '{name}' is not in subfamilies.tsv")
        if not os.path.isfile(os.path.join(config["OUT_DIR"], "subfamilies", name, "tips.tsv")):
            fail(f"no tip list for subfamily '{name}'; run scripts/build_subfamilies.py {args.family}")
out_dir = config["OUT_DIR"]
tree_file, tree_method = TREES[args.tree]
tree_dir = os.path.join(out_dir, "subfamilies", args.subfamily) if args.subfamily else os.path.join(out_dir, "tree")
out_stem = f"{args.family}_{args.subfamily}" if args.subfamily else args.family
out_path = args.out or os.path.join(out_dir, "figures", f"{out_stem}_{args.tree}.svg")
tree_path = os.path.join(tree_dir, tree_file)
title = f"{config['NAME']}, subfamily {args.subfamily}" if args.subfamily else config["NAME"]
if not os.path.isfile(tree_path):
    fail(f"'{tree_path}' not found")

tips_table = {t["tip_id"]: t for t in read_tsv(os.path.join(tree_dir, "tips.tsv"))}
references = {r["reference_id"]: r for r in read_tsv(os.path.join(config["FAMILY_DIR"], "references.tsv"))}
outgroups = read_tsv(os.path.join(config["FAMILY_DIR"], "outgroups.tsv"), required=False)
outgroup_by_id = {o["outgroup_id"]: o for o in outgroups}
tree = read_tree(tree_path)
if set(tree.tips) != set(tips_table):
    fail(f"the tips of {tree_path} differ from {tree_dir}/tips.tsv; rebuild the tree or the table")

# 1. Root
if args.subfamily:
    root_tips = [t for t, info in tips_table.items() if info["set"] != "ingroup"]
else:
    root_tips = [o["outgroup_id"] for o in outgroups if o["root"] == "TRUE"]
if root_tips and all(t in tree.tips for t in root_tips) and tree.root_on_clade(root_tips):
    names = [re.sub(r"^[^-]*-", "", t) for t in root_tips]
    root_text = "rooted on its outgroups" if args.subfamily else f"rooted on {' + '.join(names)}"
    root_key = ",".join(root_tips)
else:
    if root_tips:
        print(f"  note: the {len(root_tips)} root outgroups do not form a clade; using the midpoint root")
    tree = read_tree(tree_path)
    tree.midpoint_root()
    root_text = "midpoint root"
    root_key = "midpoint"

# 2. Scores and coverage per tip, for each HMM
tip_refs = {t: [r for r in info["references"].split(",") if r] for t, info in tips_table.items()}
tracks = []
for hmm in config["HMMS"]:
    hits = {}
    for target in ("cohort", "references"):
        path = os.path.join(out_dir, "search", f"{hmm.lower()}_{target}.domtbl")
        if not os.path.isfile(path):
            fail(f"'{path}' not found; run scripts/search_family.sbatch {args.family}")
        hits[target] = read_domtbl(path)
    per_tip = {}
    for tip in tips_table:
        hit = hits["cohort"].get(tip)
        if not hit:  # orphans carry a reference's own sequence
            found = [hits["references"][r] for r in tip_refs[tip] if r in hits["references"]]
            hit = max(found, key=lambda h: h["score"]) if found else None
        if hit:
            per_tip[tip] = hit
    if not per_tip:
        fail(f"no tip of the tree has a hit to {hmm}")
    scores = sorted(h["score"] for h in per_tip.values())
    top_score = scores[min(int(SCORE_QUANTILE * len(scores)), len(scores) - 1)]
    tracks.append(dict(hmm=hmm, hits=per_tip, top_score=top_score))

# 3. Geometry: tips evenly around the circle, the rings just outside nearly all of them
tip_order = [n.name for n in tree.preorder() if n.name]
n_tips = len(tip_order)
tip_index = {t: i for i, t in enumerate(tip_order)}
depth = tree.depths()
tip_depths = sorted(depth[id(tree.tips[t])] for t in tip_order)
reach = tip_depths[min(int(TIP_QUANTILE * n_tips), n_tips - 1)] * 1.02


def tip_angle(position):
    return 2 * math.pi * (position + 0.5) / n_tips


def polar(radius, angle):
    return WIDTH / 2 + radius * math.cos(angle), WIDTH / 2 + radius * math.sin(angle)


def xy(radius, angle):
    x, y = polar(radius, angle)
    return f"{x:.2f} {y:.2f}"


# 4. Names, and how far in the rings must stop to leave them room
marked = set(config.get("MARK_PORTALS", "").split())
labels = []
for tip, info in tips_table.items():
    in_outgroup = info.get("role") == "outgroup"
    if info["kind"] == "outgroup":
        organism = outgroup_by_id[tip]["organism"] if tip in outgroup_by_id else tip.split("-")[0]
        labels.append((tip, f"{species_prefix(organism)} {tip.split('-', 1)[1]}", True))
        continue
    refs = [references[r] for r in tip_refs[tip] if r in references and references[r]["evidence"] in LABEL_EVIDENCE]
    if refs:
        names = []
        for ref in refs:
            name = re.sub(r"\s*\(.*?\)", "", ref["protein_name"]).strip()
            if name.lower() not in [x.lower() for x in names]:
                names.append(name)
        labels.append((tip, f"{species_prefix(refs[0]['species'])} {'/'.join(names)}", in_outgroup))
    elif in_outgroup:
        labels.append((tip, tip, True))
labels.sort(key=lambda x: tip_index[x[0]])

font = FONT_SIZE
char_width = 0.5
longest = max((len(text) for _, text, _ in labels), default=0)
label_radius = OUTER_RADIUS - longest * char_width * font
if labels:
    fitting = 2 * math.pi * label_radius * 0.97 / (len(labels) * LINE_HEIGHT)
    font = min(font, fitting)
    label_radius = OUTER_RADIUS - longest * char_width * font
    if font < MIN_FONT_SIZE:
        print(f"  WARNING: {len(labels)} names only fit at {font:.1f} pt; they will be hard to read")
gap = font * LINE_HEIGHT / label_radius * n_tips / (2 * math.pi)  # the space a name needs, in tips
if labels:
    placed = spread([tip_index[t] for t, _, _ in labels], gap, 0, n_tips - 1)
else:
    placed = []

n_rings = 1 + 2 * len(tracks)
ring_width = (BAND - RING_GAP * (n_rings - 1)) / n_rings
tick_width = TICK_TRACK if marked else 0.0
ring_start = label_radius - LEADER - tick_width - BAND
scale = ring_start / reach
long_branches = sum(d > reach for d in tip_depths)

svg = []


def sector(r1, r2, start, end):
    """Closed path of a ring segment from one tip position to another."""
    if end - start + 1 == n_tips:  # SVG cannot draw a full circle as one arc
        middle = n_tips // 2
        return sector(r1, r2, start, middle - 1) + sector(r1, r2, middle, end)
    a, b = tip_angle(start - 0.5), tip_angle(end + 0.5)
    large = 1 if b - a > math.pi else 0
    return (f"M{xy(r2, a)}A{r2:.2f} {r2:.2f} 0 {large} 1 {xy(r2, b)}L{xy(r1, b)}"
            f"A{r1:.2f} {r1:.2f} 0 {large} 0 {xy(r1, a)}Z")


def draw_ring(radius, colour_of):
    """One ring: runs of tips of one colour are merged, and each colour is one path."""
    paths = {}
    start = 0
    for i in range(1, n_tips + 1):
        if i == n_tips or colour_of(tip_order[i]) != colour_of(tip_order[start]):
            colour = colour_of(tip_order[start])
            if colour:
                paths.setdefault(colour, []).append(sector(radius, radius + ring_width, start, i - 1))
            start = i
    for colour, parts in sorted(paths.items()):
        svg.append(f'<path fill="{colour}" stroke="{colour}" stroke-width="0.05" d="{"".join(parts)}"/>')


# 5. The tree, as one path; branches beyond the rings are cut at them
angle_of = {}
for node in reversed(tree.preorder()):
    if node.name:
        angle_of[id(node)] = tip_angle(tip_index[node.name])
    else:
        angles = [angle_of[id(c)] for c in node.children]
        angle_of[id(node)] = (angles[0] + angles[-1]) / 2
branches = []
limit = ring_start - 0.5
for node in tree.preorder():
    r = min(depth[id(node)] * scale, limit)
    if node.parent:
        branches.append(f"M{xy(min(depth[id(node.parent)] * scale, limit), angle_of[id(node)])}L{xy(r, angle_of[id(node)])}")
    if node.children and r > 0:
        a, b = angle_of[id(node.children[0])], angle_of[id(node.children[-1])]
        large = 1 if b - a > math.pi else 0
        branches.append(f"M{xy(r, a)}A{r:.2f} {r:.2f} 0 {large} 1 {xy(r, b)}")
root_dot = f'<circle cx="{WIDTH / 2}" cy="{WIDTH / 2}" r="3" fill="{ROOT_COLOUR}"/>'
branch_width = min(0.5, max(0.12, 300 / n_tips))  # thin lines for a whole family, heavier for a subfamily
tree_parts = [f'<path fill="none" stroke="#404040" stroke-width="{branch_width:.2f}" d="{"".join(branches)}"/>', root_dot]

# The node figure's tree: only the numbered clades, by size; a clade's radius comes from how many numbered
# clades nest inside it, so that small clades sit outward, where there is room
first_tip, last_tip = {}, {}
for node in reversed(tree.preorder()):
    if node.name:
        first_tip[id(node)] = last_tip[id(node)] = tip_index[node.name]
    else:
        first_tip[id(node)] = min(first_tip[id(c)] for c in node.children)
        last_tip[id(node)] = max(last_tip[id(c)] for c in node.children)
numbered = [n for n in tree.preorder() if n.children and n is not tree.root
            and last_tip[id(n)] - first_tip[id(n)] + 1 >= args.min_node_tips]
numbered.sort(key=lambda n: (first_tip[id(n)] - last_tip[id(n)], first_tip[id(n)]))
number_of = {id(n): i for i, n in enumerate(numbered, 1)}
levels = {}
for node in reversed(tree.preorder()):
    below = max((levels[id(c)] for c in node.children), default=0)
    levels[id(node)] = below + 1 if id(node) in number_of else below
top_level = max((levels[id(n)] for n in numbered), default=1)
step = (ring_start - 30) / top_level


def node_radius(node):
    return ring_start - 12 - (levels[id(node)] - 1) * step


def node_angle(node):
    return tip_angle((first_tip[id(node)] + last_tip[id(node)]) / 2)


def numbered_ancestor(node):
    node = node.parent
    while node is not tree.root and id(node) not in number_of:
        node = node.parent
    return node


support_labels, support_low, support_high, read_support = SUPPORT_CLASSES[args.tree]


def support_colour(node):
    try:
        value = read_support(node.label)
    except ValueError:
        return "#ffffff"
    return SUPPORT_COLOURS[0 if value < support_low else 1 if value < support_high else 2]


children_of = {}
for node in numbered:
    children_of.setdefault(id(numbered_ancestor(node)), []).append(node)
wedges, lines, circles = [], [], []
for node in numbered:
    r, angle = node_radius(node), node_angle(node)
    above = numbered_ancestor(node)
    lines.append(f"M{xy(0 if above is tree.root else node_radius(above), angle)}L{xy(r, angle)}")
    kids = children_of.get(id(node))
    if kids:
        ends = sorted([angle] + [node_angle(k) for k in kids])
        if ends[-1] > ends[0]:
            large = 1 if ends[-1] - ends[0] > math.pi else 0
            lines.append(f"M{xy(r, ends[0])}A{r:.2f} {r:.2f} 0 {large} 1 {xy(r, ends[-1])}")
    else:  # a clade with no numbered clade inside: a wedge out to the rings shows which tips it holds
        wedges.append(sector(r, ring_start - 0.5, first_tip[id(node)], last_tip[id(node)]))
    x, y = polar(r, angle)
    circles.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{NODE_RADIUS}" fill="{support_colour(node)}" '
                   f'stroke="#404040" stroke-width="0.2"/>')
    circles.append(f'<text x="{x:.2f}" y="{y:.2f}" font-size="{NODE_FONT_SIZE}" text-anchor="middle" '
                   f'dominant-baseline="central">{number_of[id(node)]}</text>')
node_parts = [f'<path fill="#ebebeb" d="{"".join(wedges)}"/>',
              f'<path fill="none" stroke="#404040" stroke-width="0.3" d="{"".join(lines)}"/>'] + circles + [root_dot]

# 6. The rings
og_colours = orthogroup_colours(list(tips_table.values()))


def orthogroup_colour(tip):
    og = tips_table[tip]["orthogroup"]
    return og_colours.get(og, OTHER_ORTHOGROUP) if og else ""


radius = ring_start
draw_ring(radius, orthogroup_colour)
for track in tracks:
    hits, top_score = track["hits"], track["top_score"]

    def binned(x, ramp):
        return ramp_colour(ramp, round(x * (HEAT_BINS - 1)) / (HEAT_BINS - 1))

    radius += ring_width + RING_GAP
    draw_ring(radius, lambda t: binned(hits[t]["score"] / top_score, SCORE_RAMP) if t in hits else NO_SCORE)
    radius += ring_width + RING_GAP
    draw_ring(radius, lambda t: binned(hits[t]["coverage"], COVERAGE_RAMP) if t in hits else NO_SCORE)
ring_end = radius + ring_width

marked_tips = [t for t in tip_order if tips_table[t]["portal"] in marked]
if marked_tips:
    ticks = "".join(f"M{xy(ring_end + 1, tip_angle(tip_index[t]))}L{xy(ring_end + tick_width - 1, tip_angle(tip_index[t]))}"
                    for t in marked_tips)
    svg.append(f'<path fill="none" stroke="#000000" stroke-width="0.5" d="{ticks}"/>')

# 7. Names and leader lines; a leader's radius grows along its way, so it never crosses the rings
leader_start = ring_end + tick_width + 1
leader_mid = leader_start + LEADER * 0.25
leader_end = label_radius - 1.5
leaders = []
texts = []
for (tip, text, is_outgroup), position in zip(labels, placed):
    a, b = tip_angle(tip_index[tip]), tip_angle(position)
    steps = max(1, int(abs(b - a) / math.radians(1.5)))
    points = [xy(leader_start, a), xy(leader_mid, a)]
    points += [xy(leader_mid + (leader_end - leader_mid) * (i / steps), a + (b - a) * (i / steps)) for i in range(1, steps + 1)]
    leaders.append("M" + "L".join(points))
    x, y = polar(label_radius, b)
    degrees = math.degrees(b)
    flip = 90 < degrees < 270
    rotation = degrees + 180 if flip else degrees
    colour = ROOT_COLOUR if is_outgroup else "#000000"
    anchor = "end" if flip else "start"
    texts.append(f'<text x="{x:.2f}" y="{y:.2f}" font-size="{font:.2f}" fill="{colour}" text-anchor="{anchor}" '
                 f'dominant-baseline="central" transform="rotate({rotation:.2f} {x:.2f} {y:.2f})">{escape(text)}</text>')
if leaders:
    svg.append(f'<path fill="none" stroke="#808080" stroke-width="0.15" d="{"".join(leaders)}"/>')
svg.extend(texts)

# 8. Legend and caption
y = WIDTH + 12
legend = []
counts = {}
for info in tips_table.values():
    counts[info["orthogroup"]] = counts.get(info["orthogroup"], 0) + 1


def swatch(x, y, colour, text):
    legend.append(f'<rect x="{x:.1f}" y="{y - 7:.1f}" width="9" height="9" fill="{colour}"/>')
    legend.append(f'<text x="{x + 12:.1f}" y="{y:.1f}" font-size="8">{escape(text)}</text>')


def colour_bar(x, y, ramp, title, low, high, mark=None):
    legend.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="8">{escape(title)}</text>')
    bar_width = 150.0
    for i in range(HEAT_BINS):
        colour = ramp_colour(ramp, i / (HEAT_BINS - 1))
        legend.append(f'<rect x="{x + bar_width * i / HEAT_BINS:.2f}" y="{y + 3:.1f}" width="{bar_width / HEAT_BINS + 0.2:.2f}" '
                      f'height="7" fill="{colour}"/>')
    legend.append(f'<text x="{x:.1f}" y="{y + 18:.1f}" font-size="7">{escape(low)}</text>')
    legend.append(f'<text x="{x + bar_width:.1f}" y="{y + 18:.1f}" font-size="7" text-anchor="end">{escape(high)}</text>')
    if mark is not None:
        mark_x = x + bar_width * mark[0]
        legend.append(f'<path d="M{mark_x:.1f} {y + 2.5:.1f}L{mark_x:.1f} {y + 11.5:.1f}" stroke="#000000" stroke-width="0.8"/>')
        legend.append(f'<text x="{mark_x:.1f}" y="{y + 18:.1f}" font-size="7" text-anchor="middle">{escape(mark[1])}</text>')


legend.append(f'<text x="30" y="{y:.1f}" font-size="9" font-weight="bold">Orthogroups (tips)</text>')
entries = [(colour, f"{og} ({counts[og]})") for og, colour in og_colours.items()]
others = sum(n for og, n in counts.items() if og and og not in og_colours)
if others:
    entries.append((OTHER_ORTHOGROUP, f"other orthogroups ({others})"))
per_column = 8
for i, (colour, text) in enumerate(entries):
    swatch(30 + 125 * (i // per_column), y + 13 + 12 * (i % per_column), colour, text)
legend_bottom = y + 13 + 12 * min(per_column, len(entries))

x_right = 30 + 125 * ((max(len(entries), 1) - 1) // per_column + 1) + 15
row = y
legend.append(f'<text x="{x_right}" y="{row:.1f}" font-size="9" font-weight="bold">HMM rings, from the centre outward</text>')
row += 14
threshold = config["MIN_HMM_COVERAGE"]
for track in tracks:
    colour_bar(x_right, row, SCORE_RAMP, f"{track['hmm']} full-sequence bit score", "0", f"≥ {track['top_score']:.0f}")
    colour_bar(x_right + 170, row, COVERAGE_RAMP, f"{track['hmm']} coverage of the HMM", "0%", "100%",
               mark=(threshold, f"{threshold:.0%}"))
    row += 30
swatch(x_right, row + 4, NO_SCORE, "no hit to this HMM (outgroups are not searched)")
row += 14
if marked_tips:
    tax = taxonomy()
    who = ", ".join(sorted(" ".join(tax[p]["species"].split()[:2]) if p in tax else p for p in marked))
    legend.append(f'<path d="M{x_right + 4} {row - 7:.1f}L{x_right + 4} {row + 2:.1f}" stroke="#000000" stroke-width="1"/>')
    legend.append(f'<text x="{x_right + 12}" y="{row:.1f}" font-size="8">{escape(who)} proteins ({len(marked_tips)})</text>')
    row += 14
legend_bottom = max(legend_bottom, row)



def write_figure(path, parts, extra_legend, caption, bottom):
    for i, line in enumerate(caption):
        extra_legend.append(f'<text x="{WIDTH / 2}" y="{bottom + 8 + 10 * i:.1f}" font-size="7.5" text-anchor="middle">'
                            f'{escape(line)}</text>')
    height = bottom + 8 + 10 * len(caption) + 6
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f'<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}pt" '
                f'height="{height:.0f}pt" viewBox="0 0 {WIDTH} {height:.0f}" font-family="Arial, Helvetica, sans-serif">\n')
        f.write(f'<rect width="{WIDTH}" height="{height:.0f}" fill="#ffffff"/>\n')
        f.write("\n".join(parts + svg + legend + extra_legend) + "\n</svg>\n")


write_figure(out_path, tree_parts, [],
             [f"{tree_method} tree of {title}, {n_tips:,} tips, {root_text} (red dot)."
              f" Names: characterized or genetically supported references (black) and outgroups (red).",
              f"Each protein's score and coverage add up the pieces of a split domain. "
              + (f"{long_branches} long terminal branches are cut at the rings." if long_branches else "")],
             legend_bottom)

# 9. The node figure and the node table
if not args.subfamily:
    node_legend = [f'<text x="{x_right}" y="{row + 2:.1f}" font-size="8">Node support:</text>']
    for i, (colour, text) in enumerate(zip(SUPPORT_COLOURS, support_labels)):
        x = x_right + 68 + 75 * i
        node_legend.append(f'<circle cx="{x}" cy="{row - 1:.1f}" r="4" fill="{colour}" stroke="#404040" stroke-width="0.2"/>')
        node_legend.append(f'<text x="{x + 7}" y="{row + 2:.1f}" font-size="8">{escape(text)}</text>')
    stem, extension = os.path.splitext(out_path)
    nodes_figure = f"{stem}_nodes{extension}"
    write_figure(nodes_figure, node_parts, node_legend,
                 [f"{tree_method} tree of {config['NAME']}, {root_text}: the {len(numbered)} clades of at least "
                  f"{args.min_node_tips} tips, numbered by size. Positions show how clades nest, not branch length.",
                  "Grey wedges show the tips of clades with no numbered clade inside. "
                  "Support is a flag, not a filter: low support means looking harder at that clade's own tree."],
                 max(legend_bottom, row + 10))

    with open(tree_path, "rb") as f:
        tree_md5 = hashlib.md5(f.read()).hexdigest()
    name_of = {tip: text for tip, text, _ in labels}
    table = []
    for node in numbered:
        inside = tree.tips_under(node)
        groups = {}
        for tip in inside:
            og = tips_table[tip]["orthogroup"] or "(none)"
            groups[og] = groups.get(og, 0) + 1
        main = sorted(groups, key=lambda og: (-groups[og], og))[:3]
        table.append([number_of[id(node)], len(inside), node.label, ", ".join(f"{og} {groups[og]}" for og in main),
                      round(groups[main[0]] / len(inside), 3), ", ".join(name_of[t] for t in inside if t in name_of),
                      sum(tips_table[t]["portal"] in marked for t in inside), inside[0], inside[-1], tree_md5, root_key])
    nodes_table = os.path.join(out_dir, "figures", f"{args.family}_{args.tree}_nodes.tsv")
    write_tsv(nodes_table, ["node", "tips", "support", "main_orthogroups", "main_share", "named", "marked", "tip_a", "tip_b",
                            "tree_md5", "root"], table)
    print(f"{nodes_figure} and {nodes_table}: {len(numbered)} clades of at least {args.min_node_tips} tips")
# 10. The figure with the subfamilies as slices
if args.subfamilies:
    slice_parts, slice_labels = [], []
    covered = {}
    for name in chosen:
        # the colour follows the subfamily's row, so it stays the same whichever subfamilies are drawn
        slot = planned.index(name)
        colour = SLICE_COLOURS[slot] if slot < len(SLICE_COLOURS) else SLICE_EXTRA_COLOUR
        ingroup = [t["tip_id"] for t in read_tsv(os.path.join(out_dir, "subfamilies", name, "tips.tsv")) if t["set"] == "ingroup"]
        positions = sorted(tip_index[t] for t in ingroup if t in tip_index)
        if len(positions) < len(ingroup):
            print(f"  WARNING: {name}: {len(ingroup) - len(positions)} ingroup tips are not in this tree")
        if not positions:
            continue
        for position in positions:
            covered[position] = covered.get(position, 0) + 1
        pieces = [[positions[0], positions[0]]]
        for position in positions[1:]:
            if position == pieces[-1][1] + 1:
                pieces[-1][1] = position
            else:
                pieces.append([position, position])
        wedges, arcs = [], []
        for first, last in pieces:
            start = min(depth[id(tree.mrca(tip_order[first:last + 1]))] * scale, limit - 2)
            wedges.append(sector(start, limit, first, last))
            arcs.append(sector(leader_start - 0.5, leader_start - 0.5 + SLICE_ARC, first, last))
        slice_parts.append(f'<path fill="{colour}" fill-opacity="{SLICE_OPACITY}" d="{"".join(wedges)}"/>')
        slice_parts.append(f'<path fill="{colour}" d="{"".join(arcs)}"/>')
        first, last = max(pieces, key=lambda piece: piece[1] - piece[0])
        angle = tip_angle((first + last) / 2)
        x, y = polar(limit - 3, angle)
        degrees = math.degrees(angle)
        flip = 90 < degrees < 270
        rotation = degrees + 180 if flip else degrees
        box = len(name) * 0.56 * SLICE_FONT_SIZE
        transform = f'transform="rotate({rotation:.2f} {x:.2f} {y:.2f})"'
        slice_labels.append(f'<rect x="{x if flip else x - box:.2f}" y="{y - SLICE_FONT_SIZE * 0.65:.2f}" width="{box:.2f}" '
                            f'height="{SLICE_FONT_SIZE * 1.3:.2f}" fill="#ffffff" fill-opacity="0.8" {transform}/>')
        slice_labels.append(f'<text x="{x:.2f}" y="{y:.2f}" font-size="{SLICE_FONT_SIZE}" text-anchor="{"start" if flip else "end"}" '
                            f'dominant-baseline="central" {transform}>{escape(name)}</text>')
        sizes = sorted((last - first + 1 for first, last in pieces), reverse=True)
        print(f"  {name}: {len(positions)} tips in {len(pieces)} piece(s)" + (f": {', '.join(map(str, sizes))}" if len(sizes) > 1 else ""))
    shared = sum(count > 1 for count in covered.values())
    outside = n_tips - len(covered)
    print(f"  {shared} tips are in more than one subfamily; {outside} tips are in none of those drawn")
    stem, extension = os.path.splitext(out_path)
    slices_figure = f"{stem}_subfamilies{extension}"
    write_figure(slices_figure, slice_parts + tree_parts + slice_labels, [],
                 [f"{tree_method} tree of {title}, {n_tips:,} tips, {root_text} (red dot). Shaded slices: subfamily "
                  f"ingroups, from their common ancestor; {outside:,} tips are in none of them.",
                  "Names: characterized or genetically supported references (black) and outgroups (red). "
                  "A subfamily in several pieces is not a clade in this tree; darker areas are in two subfamilies."],
                 legend_bottom)
    print(f"{slices_figure}: {len(chosen)} subfamilies")
print(f"{out_path}: {n_tips} tips, {len(labels)} names at {font:.1f} pt, {root_text}, "
      f"{len(marked_tips)} marked proteins, {long_branches} long branches cut")
