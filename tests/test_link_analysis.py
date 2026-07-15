"""Tests for link_analysis.py -- the Cytoscape.js element/style builder for
the force-directed Regulators + Genes view."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from network_builder import build_network  # noqa: E402
from link_analysis import (  # noqa: E402
    MULTI_GENE_LABEL,
    build_elements,
    build_styles,
    _gene_group_label,
)
from viz import REGULATOR_COLORS  # noqa: E402

DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "Connection_Matrix.xlsx")


@pytest.fixture(scope="module")
def net():
    return build_network(DATA_PATH)


def test_one_node_per_regulator_and_gene(net):
    elements = build_elements(net)
    n_regulators = len(net.regulator_labels)
    n_genes = net.T.shape[1]
    assert len(elements["nodes"]) == n_regulators + n_genes


def test_one_edge_per_regulator_gene_targeting_relationship(net):
    elements = build_elements(net)
    expected = int(net.T.values.sum())  # total number of 1s in the binary T matrix
    assert len(elements["edges"]) == expected


def test_every_edge_references_an_existing_node(net):
    elements = build_elements(net)
    node_ids = {n["data"]["id"] for n in elements["nodes"]}
    for edge in elements["edges"]:
        assert edge["data"]["source"] in node_ids
        assert edge["data"]["target"] in node_ids


def test_single_regulator_genes_grouped_with_their_regulator(net):
    elements = build_elements(net)
    singly_regulated = [g for g in net.T.columns if net.T[g].sum() == 1][0]
    reg = net.T.index[net.T[singly_regulated].astype(bool)][0]
    gene_node = next(n for n in elements["nodes"] if n["data"]["id"] == f"gene::{singly_regulated}")
    assert gene_node["data"]["label"] == _gene_group_label(reg)


def test_dual_regulator_genes_use_multi_label_and_two_edges(net):
    elements = build_elements(net)
    dual = [g for g in net.T.columns if net.T[g].sum() == 2]
    assert dual, "expected at least one dual-regulated gene in this export"
    gene = dual[0]
    gene_node = next(n for n in elements["nodes"] if n["data"]["id"] == f"gene::{gene}")
    assert gene_node["data"]["label"] == MULTI_GENE_LABEL
    matching_edges = [e for e in elements["edges"] if e["data"]["target"] == f"gene::{gene}"]
    assert len(matching_edges) == 2


def test_styles_cover_every_regulator_with_a_distinct_color(net):
    node_styles, edge_styles = build_styles(net)
    node_colors = {s.dump()["style"]["background-color"] for s in node_styles if "MULTI" not in s.dump()["selector"]}
    assert node_colors == set(REGULATOR_COLORS.values())
    # one regulator-node style + one gene-group style per regulator, plus the MULTI gene style
    assert len(node_styles) == 2 * len(net.regulator_labels) + 1
    assert len(edge_styles) == len(net.regulator_labels)


def test_regulator_node_style_bigger_than_gene_style(net):
    node_styles, _ = build_styles(net)
    by_selector = {s.dump()["selector"]: s.dump()["style"] for s in node_styles}
    reg_style = by_selector["node[label='Wcc']"]
    gene_style = by_selector[f"node[label='{_gene_group_label('Wcc')}']"]
    assert reg_style["width"] > gene_style["width"]
    assert reg_style["background-color"] == gene_style["background-color"]
