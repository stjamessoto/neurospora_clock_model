"""Tests for viz.py, focused on the 3D hub network figure -- Plotly's
Scatter3d has stricter constraints than 2D Scatter (marker.opacity and
marker.line.width must be scalars, not per-point arrays), which is an easy
thing to violate by accident when adapting 2D plotting code."""
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from network_builder import build_network  # noqa: E402
from viz import hub_network_figure, hub_network_figure_3d  # noqa: E402

DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "Connection_Matrix.xlsx")
ENRICHMENT_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "processed", "enrichment_results.csv")


@pytest.fixture(scope="module")
def net():
    return build_network(DATA_PATH)


@pytest.fixture(scope="module")
def enrichment_df():
    if not os.path.exists(ENRICHMENT_PATH):
        pytest.skip("enrichment_results.csv not cached; run src/enrichment.py first")
    return pd.read_csv(ENRICHMENT_PATH)


def test_3d_figure_builds_with_no_focus(net, enrichment_df):
    fig = hub_network_figure_3d(net, enrichment_df)
    fig.to_dict()  # triggers Plotly's schema validation
    assert len(fig.data) > 0


def test_3d_figure_builds_with_regulator_focus(net, enrichment_df):
    fig = hub_network_figure_3d(net, enrichment_df, focus_label="PostReg5 NCU08295")
    fig.to_dict()


def test_3d_figure_builds_with_wcc_focus(net, enrichment_df):
    """WCC is both a node and the only member of the 'focus' node group when
    focused on itself -- the wcc_idxs group must come out empty rather than
    erroring, since it excludes the focused node by construction."""
    fig = hub_network_figure_3d(net, enrichment_df, focus_label="Wcc")
    fig.to_dict()


def test_3d_figure_builds_without_enrichment_data(net):
    fig = hub_network_figure_3d(net, None)
    fig.to_dict()
    # spokes (10) + node trace(s) + 3 legend dummies, no chord traces at all
    assert len(fig.data) == 10 + 1 + 3


def test_3d_figure_scalar_only_marker_fields(net, enrichment_df):
    """Regression guard: Scatter3d rejects list-valued marker.opacity and
    marker.line.width (unlike 2D Scatter) -- assert every marker trace uses
    scalars for both, so a future edit can't silently reintroduce arrays."""
    fig = hub_network_figure_3d(net, enrichment_df, focus_label="PostReg5 NCU08295")
    for trace in fig.data:
        if trace.type != "scatter3d" or trace.marker is None:
            continue
        if trace.marker.opacity is not None:
            assert isinstance(trace.marker.opacity, (int, float))
        if trace.marker.line is not None and trace.marker.line.width is not None:
            assert isinstance(trace.marker.line.width, (int, float))


def test_3d_and_2d_figures_have_same_node_count(net, enrichment_df):
    fig_2d = hub_network_figure(net, enrichment_df)
    fig_3d = hub_network_figure_3d(net, enrichment_df)

    def total_markers(fig):
        return sum(len(t.x) for t in fig.data if t.mode and "markers" in t.mode and t.x and t.x[0] is not None)

    assert total_markers(fig_2d) == total_markers(fig_3d) == len(net.regulator_labels)


def test_3d_figure_class_rings_separate_on_z_axis(net):
    """The whole point of the 3D layout: transcriptional and post-
    transcriptional regulators should sit at different z heights, and WCC
    at the origin."""
    fig = hub_network_figure_3d(net, None)
    node_trace = next(t for t in fig.data if t.mode and "markers+text" == t.mode)
    zs = list(node_trace.z)
    assert 0.0 in zs  # WCC at the origin
    assert any(z > 0 for z in zs)
    assert any(z < 0 for z in zs)
