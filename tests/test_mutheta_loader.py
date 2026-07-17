"""Tests for mutheta_loader.py -- the loader for the user-supplied,
CR-delimited MuThetaDataSim_m4.txt simulated dataset."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from network_builder import REGULATOR_ORDER  # noqa: E402
from mutheta_loader import (  # noqa: E402
    N_MU,
    N_THETA,
    POST_TRANSCRIPTIONAL_LABELS,
    load_mutheta_sim,
)


@pytest.fixture(scope="module")
def df():
    return load_mutheta_sim()


def test_record_count(df):
    # Splitting the raw CR-delimited file on index resets yields 2418
    # records -- fixed by the source file, not derived from other data.
    assert len(df) == 2418


def test_columns(df):
    expected = [f"mu_{i}" for i in range(N_MU)] + [f"theta_{i}" for i in range(N_THETA)] + ["dominant_theta"]
    assert list(df.columns) == expected


def test_no_missing_values(df):
    assert not df.isna().any().any()


def test_theta_block_sums_to_one(df):
    theta_cols = [f"theta_{i}" for i in range(N_THETA)]
    sums = df[theta_cols].sum(axis=1)
    assert sums.round(6).eq(1.0).all()


def test_theta_values_are_valid_probabilities(df):
    theta_cols = [f"theta_{i}" for i in range(N_THETA)]
    assert (df[theta_cols].values >= 0).all()
    assert (df[theta_cols].values <= 1).all()


def test_mu_block_roughly_0_to_100(df):
    # Not a hard biological/statistical bound, just documents the observed
    # range so a future format change (e.g. a units switch) trips a test.
    mu_cols = [f"mu_{i}" for i in range(N_MU)]
    assert (df[mu_cols].values >= 0).all()
    assert (df[mu_cols].values <= 100).all()


def test_dominant_theta_is_one_of_the_theta_columns(df):
    theta_cols = set(f"theta_{i}" for i in range(N_THETA))
    assert set(df["dominant_theta"].unique()) <= theta_cols


def test_theta_0_dominates_far_more_than_uniform(df):
    """Structural fact that motivates (but does not prove) the resemblance
    to this project's real post-transcriptional regulator set: theta_0 wins
    the argmax much more often than a uniform 1/6, mirroring lhp-1's
    outsized real target share among the 6 post-transcriptional regulators
    (see mutheta_loader module docstring)."""
    share = df["dominant_theta"].value_counts(normalize=True)
    assert share["theta_0"] > 2 * (1 / N_THETA)


def test_post_transcriptional_labels_match_network_builder():
    assert POST_TRANSCRIPTIONAL_LABELS == REGULATOR_ORDER[-6:]
    assert len(POST_TRANSCRIPTIONAL_LABELS) == N_THETA


def test_label_theta_as_post_transcriptional_option():
    labeled = load_mutheta_sim(label_theta_as_post_transcriptional=True)
    for label in POST_TRANSCRIPTIONAL_LABELS:
        assert label in labeled.columns
    assert set(labeled["dominant_theta"].unique()) <= set(POST_TRANSCRIPTIONAL_LABELS)
