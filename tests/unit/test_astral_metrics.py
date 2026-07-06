"""Unit tests for workflow/scripts/astral_metrics.py.

Covers internal-node numbering used to key the per-node support/quartet metrics.
(newick_to_nhx is intentionally not tested here: its output depends on ASTRAL's
exact quoted-annotation format and is better exercised by an integration test.)
"""
import pytest

pytest.importorskip("ete3")
from ete3 import Tree

import astral_metrics as am


def test_add_node_numbers_numbers_only_internal_nodes():
    tree = Tree("((A,B),C);")
    am.add_node_numbers(tree)

    internal = [n for n in tree.traverse() if not n.is_leaf()]
    leaves = [n for n in tree.traverse() if n.is_leaf()]

    # Every internal node gets a sequential number starting at 1; leaves get none.
    numbers = sorted(n.node_number for n in internal)
    assert numbers == list(range(1, len(internal) + 1))
    assert all(not hasattr(leaf, "node_number") for leaf in leaves)
