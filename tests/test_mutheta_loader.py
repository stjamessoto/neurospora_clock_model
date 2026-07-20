"""Tests for mutheta_loader.py -- the loader for the user-supplied,
CR-delimited MuThetaDataSim_m4.txt simulated dataset and its gene list
(ALLNCU.txt), plus the per-gene validation against real regulator identity."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from network_builder import REGULATOR_ORDER  # noqa: E402
from mutheta_loader import (  # noqa: E402
    N_MU,
    N_THETA,
    POST_TRANSCRIPTIONAL_LABELS,
    best_theta_to_regulator_mapping,
    load_gene_list,
    load_mutheta_sim,
    load_mutheta_sim_with_genes,
    real_post_transcriptional_regulator,
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


# --- ALLNCU.txt gene list ---------------------------------------------


def test_gene_list_length_matches_sim_record_count(df):
    assert len(load_gene_list()) == len(df)


def test_gene_list_has_expected_duplication():
    # Fixed fact about the source file: 2418 lines, 2216 distinct genes.
    genes = load_gene_list()
    assert len(genes) == 2418
    assert len(set(genes)) == 2216


def test_every_distinct_gene_is_a_real_target_gene():
    """Every gene in ALLNCU.txt is a real column in the cleaned
    Connection_Matrix.xlsx target matrix -- strong evidence this dataset
    was drawn from the same gene universe as the rest of this project."""
    from data_loader import load_clean

    targets, _, _ = load_clean("data/raw/Connection_Matrix.xlsx")
    genes = set(load_gene_list())
    assert genes <= set(targets.columns)


def test_load_mutheta_sim_with_genes_joins_by_row_order():
    joined = load_mutheta_sim_with_genes()
    assert len(joined) == 2418
    assert list(joined["gene"]) == load_gene_list()


# --- Per-gene validation against the real network -----------------------


def test_real_post_transcriptional_regulator_lookup():
    lookup = real_post_transcriptional_regulator()
    assert 800 < len(lookup) < 900  # 887 at time of writing; loose bound
    assert set(lookup.values()) <= set(POST_TRANSCRIPTIONAL_LABELS)


def test_theta_argmax_does_not_beat_majority_baseline_for_real_regulator_identity():
    """Locks in the module's central finding: brute-forcing all 720
    theta-column -> regulator-label permutations, the best achievable
    per-gene accuracy against each gene's real single post-transcriptional
    regulator is *below* the trivial "always guess the majority regulator"
    baseline. If this ever flips, it's a real finding worth updating the
    README over -- not something to just delete on failure."""
    result = best_theta_to_regulator_mapping()
    assert result["n_genes_evaluated"] > 500
    assert result["best_acc"] < result["majority_baseline_acc"]
    # sanity: the permutation search itself should average out to chance
    assert result["mean_acc_all_permutations"] == pytest.approx(1 / N_THETA, abs=0.01)
