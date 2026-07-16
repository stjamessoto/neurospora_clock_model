"""Tests for gene_layout_3d.py -- the cached 3D force-directed layout for
the Regulators + Genes graph. Full recomputation (~45s) is only exercised
once here; the cached result on disk is what the app actually reads."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from network_builder import build_network  # noqa: E402
from gene_layout_3d import DEFAULT_CACHE_PATH, load_result  # noqa: E402

DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "Connection_Matrix.xlsx")


@pytest.fixture(scope="module")
def net():
    return build_network(DATA_PATH)


@pytest.fixture(scope="module")
def cached_result():
    result = load_result()
    if result is None:
        pytest.skip(f"{DEFAULT_CACHE_PATH} not cached; run src/gene_layout_3d.py first")
    return result


def test_cached_layout_covers_every_node(net, cached_result):
    expected_ids = {f"reg::{lab}" for lab in net.regulator_labels} | {f"gene::{g}" for g in net.T.columns}
    assert set(cached_result.node_ids) == expected_ids


def test_cached_positions_are_3d_and_finite(cached_result):
    assert cached_result.positions.shape == (len(cached_result.node_ids), 3)
    assert np.isfinite(cached_result.positions).all()


def test_compute_layout_small_graph_deterministic_with_seed():
    """Full spring_layout is slow on the real graph (~45s) -- this exercises
    compute_layout itself on a small synthetic network so the test suite
    doesn't pay that cost, and checks the seed makes it reproducible."""
    from dataclasses import dataclass

    import pandas as pd

    from gene_layout_3d import compute_layout

    @dataclass
    class _FakeNet:
        T: pd.DataFrame
        regulator_labels: list

    T = pd.DataFrame(
        [[1, 0, 1], [0, 1, 1]],
        index=["Wcc", "Reg1 NCU07392"],
        columns=["NCU00001", "NCU00002", "NCU00003"],
    )
    net = _FakeNet(T=T, regulator_labels=["Wcc", "Reg1 NCU07392"])

    result_a = compute_layout(net, seed=1, iterations=10)
    result_b = compute_layout(net, seed=1, iterations=10)
    assert result_a.node_ids == result_b.node_ids
    assert np.allclose(result_a.positions, result_b.positions)
    assert result_a.positions.shape == (5, 3)  # 2 regulators + 3 genes
