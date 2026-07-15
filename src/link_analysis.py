"""Builds Cytoscape.js elements/styles for the force-directed regulator +
target-gene graph (Regulators + Genes subtab), rendered via
``st_link_analysis``. Each regulator gets its own node/gene/edge color from
``viz.REGULATOR_COLORS`` -- one hub, one color -- so the graph reads the
same way as the reference hub-and-spoke layout it's modeled on. Layout
itself (force-directed clustering) is computed client-side by Cytoscape.js's
"fcose" algorithm, not here; this module only supplies the graph data and
its styling.

``st_link_analysis``'s public ``NodeStyle``/``EdgeStyle`` classes don't
expose node size, and this graph needs regulator nodes to render much
larger than the ~3,026 gene dots hanging off them. ``_RawStyle`` is a
minimal stand-in with the same ``.dump()`` contract those classes use
internally, so a raw Cytoscape.js style dict (including ``width``/
``height``) can be passed through the same ``node_styles``/``edge_styles``
parameters without forking the package.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from network_builder import ClockNetwork
from viz import REGULATOR_COLORS

MULTI_GENE_LABEL = "MULTI"
MULTI_GENE_COLOR = "#898781"  # neutral ink-muted -- genes with 2 regulators, no single hub color

REGULATOR_NODE_SIZE = 46
GENE_NODE_SIZE = 9


@dataclass
class _RawStyle:
    _dump: dict[str, Any]

    def dump(self) -> dict[str, Any]:
        return self._dump


def _node_style(label: str, color: str, size: float, caption: str | None = None) -> _RawStyle:
    style: dict[str, Any] = {"background-color": color, "width": size, "height": size}
    if caption:
        style["label"] = f"data({caption})"
    return _RawStyle({"selector": f"node[label='{label}']", "style": style})


def _edge_style(label: str, color: str) -> _RawStyle:
    return _RawStyle({"selector": f"edge[label='{label}']", "style": {"line-color": color}})


def _gene_group_label(regulator_label: str) -> str:
    return f"{regulator_label}::gene"


def build_elements(net: ClockNetwork) -> dict[str, list[dict]]:
    """Cytoscape elements: one node per regulator, one node per target gene,
    one edge per (regulator, gene) targeting relationship.

    A gene's node ``label`` (used for style-group selection, not display
    text) is its single regulator's gene-group -- so it renders in that
    regulator's color -- or ``MULTI_GENE_LABEL`` for the small set of genes
    targeted by two regulators (this project's max per gene), which get one
    edge to each.
    """
    nodes = [{"data": {"id": f"reg::{lab}", "label": lab, "name": lab}} for lab in net.regulator_labels]
    edges = []

    for gene in net.T.columns:
        regs = list(net.T.index[net.T[gene].astype(bool)])
        gene_group = _gene_group_label(regs[0]) if len(regs) == 1 else MULTI_GENE_LABEL
        nodes.append({"data": {"id": f"gene::{gene}", "label": gene_group, "name": gene}})
        for reg in regs:
            edges.append({
                "data": {
                    "id": f"{reg}->{gene}",
                    "label": reg,
                    "source": f"reg::{reg}",
                    "target": f"gene::{gene}",
                }
            })

    return {"nodes": nodes, "edges": edges}


def build_styles(net: ClockNetwork) -> tuple[list[_RawStyle], list[_RawStyle]]:
    """Node/edge style lists for ``st_link_analysis``: one color per
    regulator, applied to that regulator's own (large, captioned) node, its
    exclusively-targeted genes (small, uncaptioned), and every edge leaving
    it; plus one neutral group for the dual-regulated genes."""
    node_styles = []
    edge_styles = []
    for lab in net.regulator_labels:
        color = REGULATOR_COLORS[lab]
        node_styles.append(_node_style(lab, color, REGULATOR_NODE_SIZE, caption="name"))
        node_styles.append(_node_style(_gene_group_label(lab), color, GENE_NODE_SIZE))
        edge_styles.append(_edge_style(lab, color))
    node_styles.append(_node_style(MULTI_GENE_LABEL, MULTI_GENE_COLOR, GENE_NODE_SIZE))
    return node_styles, edge_styles
