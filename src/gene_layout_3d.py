"""Cached 3D force-directed layout for the Regulators + Genes graph.

``networkx.spring_layout(dim=3)`` (Fruchterman-Reingold) on this project's
full regulator+target-gene graph (3,037 nodes, 3,255 edges) takes ~45s --
too slow to recompute on every Streamlit rerun, so it's computed once here
and cached to disk, the same pattern ``null_models.py`` and ``enrichment.py``
use for their own expensive one-time computations.
"""
from __future__ import annotations

import os
from dataclasses import dataclass

import networkx as nx
import numpy as np

from network_builder import ClockNetwork, build_network

DEFAULT_CACHE_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "processed", "gene_layout_3d.npz")


@dataclass
class GeneLayout3DResult:
    node_ids: list[str]      # "reg::<regulator label>" / "gene::<gene name>"
    positions: np.ndarray    # (N, 3), row i is node_ids[i]'s (x, y, z)


def build_graph(net: ClockNetwork) -> nx.Graph:
    G = nx.Graph()
    for lab in net.regulator_labels:
        G.add_node(f"reg::{lab}")
    for gene in net.T.columns:
        G.add_node(f"gene::{gene}")
        for lab in net.T.index[net.T[gene].astype(bool)]:
            G.add_edge(f"reg::{lab}", f"gene::{gene}")
    return G


def compute_layout(
    net: ClockNetwork, seed: int = 42, iterations: int = 80, k: float | None = 0.4
) -> GeneLayout3DResult:
    """``k`` is spring_layout's optimal node-node distance; networkx's
    default (``1/sqrt(n)`` ~= 0.018 for this graph's 3,037 nodes) packs
    regulator hubs and their gene clusters too tightly to see individual
    edges. A larger fixed ``k`` (regulators have no edges between each
    other, so their only separation comes from this repulsion term) spreads
    both the hubs apart from each other and the genes apart within each
    hub; ``iterations`` is raised alongside it so the layout actually
    converges at the wider spacing (~44s total at these settings, still
    only paid once since the result is cached).
    """
    G = build_graph(net)
    pos = nx.spring_layout(G, dim=3, seed=seed, iterations=iterations, k=k)
    node_ids = list(G.nodes())
    positions = np.array([pos[n] for n in node_ids])
    return GeneLayout3DResult(node_ids=node_ids, positions=positions)


def save_result(result: GeneLayout3DResult, path: str = DEFAULT_CACHE_PATH) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, node_ids=np.array(result.node_ids, dtype=object), positions=result.positions)


def load_result(path: str = DEFAULT_CACHE_PATH) -> GeneLayout3DResult | None:
    if not os.path.exists(path):
        return None
    data = np.load(path, allow_pickle=True)
    return GeneLayout3DResult(node_ids=list(data["node_ids"]), positions=data["positions"])


if __name__ == "__main__":
    net = build_network(os.path.join(os.path.dirname(__file__), "..", "data", "raw", "Connection_Matrix.xlsx"))
    result = compute_layout(net)
    save_result(result)
    print(f"Saved 3D layout for {len(result.node_ids)} nodes to {DEFAULT_CACHE_PATH}")
