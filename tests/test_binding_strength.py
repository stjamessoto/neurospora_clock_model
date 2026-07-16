"""Tests for binding_strength.py -- the loader for the user-supplied
function x regulator continuous binding-strength table."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from network_builder import REGULATOR_ORDER  # noqa: E402
from binding_strength import load_binding_strength_by_function  # noqa: E402


@pytest.fixture(scope="module")
def table():
    return load_binding_strength_by_function()


def test_shape(table):
    assert table.shape == (21, 11)


def test_columns_are_canonical_regulator_labels_in_order(table):
    assert list(table.columns) == REGULATOR_ORDER


def test_no_missing_values(table):
    assert not table.isna().any().any()


def test_values_are_nonnegative_and_bounded(table):
    assert (table.values >= 0).all()
    assert (table.values < 1).all()  # these are binding strengths, not raw target counts


def test_known_cell_values_spot_check(table):
    # transcribed directly from the source table -- exact spot checks, not tolerances
    assert table.loc["Ribosome", "Wcc"] == 0.0
    assert table.loc["RNA transport", "Wcc"] == pytest.approx(0.7103)
    assert table.loc["Mismatch repair", "PostReg10 NCU09349"] == 0.0
    assert table.loc["Purine metabolism", "Reg4 NCU07155"] == pytest.approx(0.7651)


def test_most_rows_have_a_true_zero(table):
    """Most function rows have at least one 0.0 cell (not the dominant
    regulator for that function) -- if this stopped being true the "blank
    means zero" heatmap encoding would be misleading. Not universal: e.g.
    "Biosynthesis of amino acids" is touched by all 11 regulators."""
    has_zero = (table == 0.0).any(axis=1)
    assert has_zero.mean() > 0.8
    assert not table.loc["Biosynthesis of amino acids"].eq(0.0).any()
