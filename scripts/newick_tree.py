#!/usr/bin/env python3
# newick_tree.py
# Description: Shared helpers for reading and querying gene family trees (not run on its own):
#              reads Newick trees with support labels (FastTree, IQ-TREE), roots an unrooted tree on
#              a tip, a clade or the midpoint, and finds clades. Written without recursion, for trees of thousands of tips.
#              A node's length and label belong to the branch above it, so they stay with the same
#              split when the tree is rooted differently.
# Author: Victor M. Gonzalez Ramos
# Date: 2026-11-05
# Version: 1.0

import sys


class Node:
    __slots__ = ("name", "length", "label", "parent", "children")

    def __init__(self, parent=None):
        self.name, self.label = "", ""
        self.length = 0.0
        self.parent = parent
        self.children = []


class Tree:
    def __init__(self, root, tips):
        self.root = root
        self.tips = tips  # tip name -> Node

    def root_on_tip(self, name):
        """Roots the tree on a tip: every clade is then defined relative to that tip."""
        path = [self.tips[name]]
        while path[-1].parent:
            path.append(path[-1].parent)
        for child, parent in zip(reversed(path[:-1]), reversed(path[1:])):
            parent.children.remove(child)
            child.children.append(parent)
            parent.parent = child
            parent.length, parent.label = child.length, child.label
        self.root = path[0]
        self.root.parent, self.root.length, self.root.label = None, 0.0, ""
        self._remove_if_pass_through(path[-1])

    def root_on_branch(self, node, distance):
        """Roots the tree on the branch above a node, at a distance above the node."""
        length, label = node.length, node.label
        above = node.parent
        above.children.remove(node)
        path = [above]
        while path[-1].parent:
            path.append(path[-1].parent)
        for child, parent in zip(reversed(path[:-1]), reversed(path[1:])):
            parent.children.remove(child)
            child.children.append(parent)
            parent.parent = child
            parent.length, parent.label = child.length, child.label
        root = Node()
        root.children = [node, above]
        node.parent, above.parent = root, root
        node.length, above.length = distance, length - distance
        node.label, above.label = label, label
        self._remove_if_pass_through(path[-1])
        self.root = root

    @staticmethod
    def _remove_if_pass_through(node):
        """The old root of a rooted tree is left with one child: it is no longer a node of the tree."""
        if len(node.children) == 1 and node.parent and not node.name:
            child = node.children[0]
            child.parent = node.parent
            child.length += node.length
            node.parent.children[node.parent.children.index(node)] = child

    def root_on_clade(self, names):
        """Roots on the branch that separates the named tips from the rest, at its middle. False, leaving
        the tree rooted on the first tip, if they do not form a clade in the unrooted tree."""
        self.root_on_tip(names[0])
        rest = [t for t in self.tips if t not in names]
        node = self.mrca(rest)
        if len(self.tips_under(node)) != len(rest):
            return False
        self.root_on_branch(node, node.length / 2)
        return True

    def midpoint_root(self):
        """Roots halfway along the longest path between two tips."""
        self.root_on_tip(self.farthest_tip())
        last = self.farthest_tip()
        half = self.depths()[id(self.tips[last])] / 2
        node, travelled = self.tips[last], 0.0
        while travelled + node.length < half:
            travelled += node.length
            node = node.parent
        self.root_on_branch(node, half - travelled)

    def preorder(self):
        """All nodes, parents before children."""
        nodes, stack = [], [self.root]
        while stack:
            current = stack.pop()
            nodes.append(current)
            stack.extend(reversed(current.children))
        return nodes

    def depths(self):
        """Distance from the root of every node, by id(node)."""
        depth = {}
        for node in self.preorder():
            depth[id(node)] = depth[id(node.parent)] + node.length if node.parent else 0.0
        return depth

    def farthest_tip(self):
        """The tip farthest from the current root."""
        depth = self.depths()
        return max(self.tips, key=lambda t: (depth[id(self.tips[t])], t))

    def mrca(self, names):
        """Most recent common ancestor of named tips (with one tip, the tip itself)."""
        ancestors = []
        node = self.tips[names[0]]
        while node:
            ancestors.append(node)
            node = node.parent
        position = {id(n): i for i, n in enumerate(ancestors)}
        lowest = 0
        for name in names[1:]:
            node = self.tips[name]
            while id(node) not in position:
                node = node.parent
            lowest = max(lowest, position[id(node)])
        return ancestors[lowest]

    def tips_under(self, node):
        """Names of the tips below a node, in tree order."""
        names, stack = [], [node]
        while stack:
            current = stack.pop()
            if current.name:
                names.append(current.name)
            stack.extend(reversed(current.children))
        return names


def read_tree(path):
    """Reads one Newick tree. Tips are nodes with a name; internal nodes keep their label (support)."""
    with open(path) as f:
        text = "".join(f.read().split())
    if not text.endswith(";") or text.count(";") != 1:
        print(f"ERROR: '{path}' should hold exactly one Newick tree", file=sys.stderr)
        sys.exit(1)
    root = node = Node()
    tips = {}
    i, end = 0, len(text) - 1
    while i < end:
        c = text[i]
        if c == "(":
            child = Node(node)
            node.children.append(child)
            node = child
            i += 1
        elif c == ",":
            child = Node(node.parent)
            node.parent.children.append(child)
            node = child
            i += 1
        elif c == ")":
            node = node.parent
            i += 1
        else:  # a name or label, then an optional branch length
            j = i
            while j < end and text[j] not in "(),":
                j += 1
            word, _, length = text[i:j].partition(":")
            if node.children:
                node.label = word
            elif word:
                node.name = word
                tips[word] = node
            if length:
                node.length = float(length)
            i = j
    return Tree(root, tips)
