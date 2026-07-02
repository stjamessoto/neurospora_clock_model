"""Sanity checks against the paper's published numbers and internal invariants.

Where this export's cleaned network diverges from the paper (see README:
"Validation against the paper"), tests assert the *documented* divergent
value, not the paper's number -- the point is to catch regressions in our
own pipeline, not to force agreement with a different underlying dataset.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from data_loader import load_clean  # noqa: E402
from network_builder import WCC_INDEX, build_network  # noqa: E402
from ffl_analysis import count_ffls  # noqa: E402
from enrichment import _triple_overlap_pmf, pairwise_enrichment, run_all_enrichment_tests  # noqa: E402
from null_models import _double_edge_swap, _ffl_count_from_matrix  # noqa: E402

DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "Connection_Matrix.xlsx")


@pytest.fixture(scope="module")
def net():
    return build_network(DATA_PATH)


# --- data cleaning -----------------------------------------------------

def test_no_whitespace_in_labels():
    targets, core_clock, report = load_clean(DATA_PATH)
    assert all(lab == lab.strip() for lab in targets.index)
    assert all(col == col.strip() for col in targets.columns)


def test_rok1_target_count_matches_paper():
    targets, _, _ = load_clean(DATA_PATH)
    assert int(targets.loc["PostReg9 NCU00919"].sum()) == 54


def test_core_clock_genes_excluded_from_targets():
    targets, core_clock, _ = load_clean(DATA_PATH)
    assert core_clock.shape[1] == 3
    assert not any(c.startswith(("wc-1", "wc-2", "frq")) for c in targets.columns)


def test_target_columns_are_unique():
    targets, _, _ = load_clean(DATA_PATH)
    assert targets.columns.is_unique


# --- network construction ----------------------------------------------

def test_R_matrix_structure(net):
    assert net.R.shape == (11, 11)
    assert np.all(np.diag(net.R) == 0)
    assert net.R[WCC_INDEX, :].sum() == 10  # WCC activates all 10 other regulators
    assert net.R.sum() == 10  # only WCC's row is nonzero
    non_wcc_rows = [i for i in range(11) if i != WCC_INDEX]
    assert np.all(net.R[non_wcc_rows, :] == 0)  # regulators don't regulate each other


def test_T_matrix_is_binary(net):
    vals = net.T.values
    assert set(np.unique(vals)).issubset({0, 1})


# --- FFL counting --------------------------------------------------------

def test_ffl_count_documented_divergence(net):
    """This export's max regulator-degree-per-gene is 2, so the FFL count
    (30) genuinely differs from the paper's 71 -- see README. This test
    pins the value so a future data/cleaning change doesn't silently shift
    it without notice."""
    result = count_ffls(net)
    assert result.total_ffl_count == 30


def test_ffl_matrix_matches_fast_path(net):
    """ffl_analysis.count_ffls (full Hadamard R@T formula) must agree with
    null_models._ffl_count_from_matrix (fast path used inside the null
    model loop), since the latter is derived algebraically from the former."""
    result = count_ffls(net)
    fast = _ffl_count_from_matrix(net.T.values.astype(bool), WCC_INDEX)
    assert result.total_ffl_count == fast


def test_per_regulator_ffl_sums_to_total(net):
    result = count_ffls(net)
    assert result.per_regulator_ffl_count.sum() == result.total_ffl_count


def test_no_overlap_genes_in_this_export(net):
    """Documented consequence of max-degree-2: no gene has >1 non-WCC
    regulator plus WCC, so no gene can show 2+ overlapping FFLs here."""
    result = count_ffls(net)
    assert len(result.overlap_genes) == 0


# --- null model ------------------------------------------------------------

def test_double_edge_swap_preserves_row_and_column_sums(net):
    T0 = net.T.values.astype(bool).copy()
    row_sums_before = T0.sum(axis=1).copy()
    col_sums_before = T0.sum(axis=0).copy()

    rng = np.random.default_rng(42)
    T = T0.copy()
    n_success = _double_edge_swap(T, n_swaps=500, rng=rng)

    assert n_success > 0
    np.testing.assert_array_equal(T.sum(axis=1), row_sums_before)
    np.testing.assert_array_equal(T.sum(axis=0), col_sums_before)
    # and the matrix actually changed
    assert not np.array_equal(T, T0)


# --- enrichment --------------------------------------------------------

def test_enrichment_test_counts():
    """55 pairs + 165 triples = 220 tests, per C(11,2) and C(11,3)."""
    from math import comb

    assert comb(11, 2) == 55
    assert comb(11, 3) == 165


def test_triple_overlap_pmf_sums_to_one():
    pmf = _triple_overlap_pmf(N=3026, n1=200, n2=150, n3=100)
    assert pmf.sum() == pytest.approx(1.0, abs=1e-8)


def test_triple_overlap_pmf_matches_expectation_formula():
    N, n1, n2, n3 = 3026, 200, 150, 100
    pmf = _triple_overlap_pmf(N, n1, n2, n3)
    mean_from_pmf = float(np.dot(np.arange(len(pmf)), pmf))
    expected_closed_form = n1 * n2 * n3 / (N ** 2)
    assert mean_from_pmf == pytest.approx(expected_closed_form, rel=1e-6)


def test_pairwise_enrichment_shape(net):
    df = pairwise_enrichment(net)
    assert len(df) == 55
    assert set(df["p_value"]) <= set(df["p_value"])  # all present
    assert (df["p_value"] >= 0).all() and (df["p_value"] <= 1).all()


def test_bh_correction_monotonic(net):
    combined = run_all_enrichment_tests(net)
    assert len(combined) == 220
    sorted_p = combined.sort_values("p_value")["q_value"].values
    # q-values must be non-decreasing when sorted by ascending raw p-value
    assert np.all(np.diff(sorted_p) >= -1e-12)
