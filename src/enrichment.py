"""Pairwise and triple-wise regulator co-targeting enrichment tests.

Novel extension beyond the paper (see README): the paper never tested
whether any two (or three) regulators share more target genes than
expected by chance. Null model: each regulator's target set is treated as
a fixed-size random subset of the N-gene universe, drawn independently of
every other regulator's target set (their sizes -- "row sums" of T -- are
taken as given, exactly as in the paper's fitted network).

- Pairwise (55 = C(11,2) tests): exact hypergeometric test (equivalent to
  a two-sided Fisher's exact test on the 2x2 contingency table for each
  regulator pair), via scipy.stats.hypergeom.
- Triple-wise (165 = C(11,3) tests): exact distribution of a 3-way
  intersection size under the same independent-random-subset null,
  computed by convolving two hypergeometric laws:
      P(|A n B| = m)              ~ Hypergeom(N, n_A, n_B)
      P(|C n (A n B)| = k | m)    ~ Hypergeom(N, m, n_C)
      P(|A n B n C| = k) = sum_m P(|A n B|=m) * P(k | m)
  This is exact (not a normal/simulation approximation) because A, B, C
  are independent draws, so the "population" being intersected with C is
  exactly the (random-sized) set A n B.

All 55 + 165 = 220 p-values (two-sided) are pooled and Benjamini-Hochberg
FDR corrected at alpha = 0.05 (statsmodels multipletests), mirroring the
significance-table format used in the companion S. cerevisiae subnetwork-
motif project.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from data_loader import REGULATOR_INFO
from network_builder import ClockNetwork, build_network

ALPHA = 0.05


def _pairwise_pvalue(observed: int, N: int, n1: int, n2: int) -> tuple[float, float]:
    """Two-sided exact hypergeometric p-value for |A n B| = observed.

    Returns (p_upper, p_two_sided) where p_upper = P(X >= observed) is the
    enrichment-direction tail.
    """
    dist = hypergeom(N, n1, n2)
    p_upper = float(dist.sf(observed - 1))
    p_lower = float(dist.cdf(observed))
    p_two_sided = min(1.0, 2 * min(p_upper, p_lower))
    return p_upper, p_two_sided


def _triple_overlap_pmf(N: int, n1: int, n2: int, n3: int) -> np.ndarray:
    """Exact pmf of |A n B n C| for independent random subsets of sizes n1,n2,n3 from N.

    Returns an array `pmf` where pmf[k] = P(|A n B n C| = k).
    """
    m_min = max(0, n1 + n2 - N)
    m_max = min(n1, n2)
    ms = np.arange(m_min, m_max + 1)
    p_m = hypergeom.pmf(ms, N, n1, n2)

    max_k = min(m_max, n3)
    pmf = np.zeros(max_k + 1)
    for m, pm in zip(ms, p_m):
        if pm <= 0:
            continue
        k_max_m = min(int(m), n3)
        ks = np.arange(0, k_max_m + 1)
        p_k_given_m = hypergeom.pmf(ks, N, int(m), n3)
        pmf[: k_max_m + 1] += pm * p_k_given_m
    return pmf


def _triple_pvalue(observed: int, N: int, n1: int, n2: int, n3: int) -> tuple[float, float, float]:
    """Two-sided exact p-value for a 3-way overlap, plus its null expectation.

    Returns (expected, p_upper, p_two_sided).
    """
    expected = n1 * n2 * n3 / (N ** 2)
    pmf = _triple_overlap_pmf(N, n1, n2, n3)
    observed = min(observed, len(pmf) - 1)  # guard, should not trigger
    p_upper = float(pmf[observed:].sum())
    p_lower = float(pmf[: observed + 1].sum())
    p_two_sided = min(1.0, 2 * min(p_upper, p_lower))
    return expected, p_upper, p_two_sided


def pairwise_enrichment(net: ClockNetwork) -> pd.DataFrame:
    T = net.T
    N = T.shape[1]
    labels = list(T.index)
    names = {lab: REGULATOR_INFO[lab]["gene"] for lab in labels}
    sizes = T.sum(axis=1)

    rows = []
    for a, b in combinations(labels, 2):
        set_a = set(T.columns[T.loc[a].values.astype(bool)])
        set_b = set(T.columns[T.loc[b].values.astype(bool)])
        observed = len(set_a & set_b)
        n1, n2 = int(sizes[a]), int(sizes[b])
        expected = n1 * n2 / N
        p_upper, p_two_sided = _pairwise_pvalue(observed, N, n1, n2)
        rows.append(
            {
                "test_type": "pair",
                "regulators": f"{names[a]} + {names[b]}",
                "reg_a": a,
                "reg_b": b,
                "reg_c": None,
                "n_a": n1,
                "n_b": n2,
                "n_c": np.nan,
                "observed_overlap": observed,
                "expected_overlap": expected,
                "direction": "enriched" if observed >= expected else "depleted",
                "p_upper": p_upper,
                "p_value": p_two_sided,
            }
        )
    return pd.DataFrame(rows)


def triple_enrichment(net: ClockNetwork) -> pd.DataFrame:
    T = net.T
    N = T.shape[1]
    labels = list(T.index)
    names = {lab: REGULATOR_INFO[lab]["gene"] for lab in labels}
    sizes = T.sum(axis=1)
    bool_sets = {lab: set(T.columns[T.loc[lab].values.astype(bool)]) for lab in labels}

    rows = []
    for a, b, c in combinations(labels, 3):
        observed = len(bool_sets[a] & bool_sets[b] & bool_sets[c])
        n1, n2, n3 = int(sizes[a]), int(sizes[b]), int(sizes[c])
        expected, p_upper, p_two_sided = _triple_pvalue(observed, N, n1, n2, n3)
        rows.append(
            {
                "test_type": "triple",
                "regulators": f"{names[a]} + {names[b]} + {names[c]}",
                "reg_a": a,
                "reg_b": b,
                "reg_c": c,
                "n_a": n1,
                "n_b": n2,
                "n_c": n3,
                "observed_overlap": observed,
                "expected_overlap": expected,
                "direction": "enriched" if observed >= expected else "depleted",
                "p_upper": p_upper,
                "p_value": p_two_sided,
            }
        )
    return pd.DataFrame(rows)


def run_all_enrichment_tests(net: ClockNetwork, alpha: float = ALPHA) -> pd.DataFrame:
    """Run all 55 pairwise + 165 triple tests, pool, and apply BH-FDR correction."""
    pairs = pairwise_enrichment(net)
    triples = triple_enrichment(net)
    combined = pd.concat([pairs, triples], ignore_index=True)

    assert len(combined) == 55 + 165, f"expected 220 tests, got {len(combined)}"

    reject, q_values, _, _ = multipletests(combined["p_value"].values, alpha=alpha, method="fdr_bh")
    combined["q_value"] = q_values
    combined["significant_bh"] = reject
    combined = combined.sort_values("q_value", ascending=True).reset_index(drop=True)
    return combined


def summarize(combined: pd.DataFrame, alpha: float = ALPHA) -> str:
    n_sig = int(combined["significant_bh"].sum())
    n_pairs = int((combined["test_type"] == "pair").sum())
    n_triples = int((combined["test_type"] == "triple").sum())
    n_sig_pairs = int(((combined["test_type"] == "pair") & combined["significant_bh"]).sum())
    n_sig_triples = int(((combined["test_type"] == "triple") & combined["significant_bh"]).sum())

    lines = [
        f"Ran {len(combined)} tests ({n_pairs} pairs + {n_triples} triples), "
        f"BH-FDR corrected at alpha={alpha}.",
        f"Significant after correction: {n_sig} total ({n_sig_pairs} pairs, {n_sig_triples} triples).",
        "",
        "Top 15 by corrected q-value:",
        combined.head(15)[
            ["regulators", "test_type", "observed_overlap", "expected_overlap", "direction", "p_value", "q_value"]
        ].to_string(index=False),
    ]
    return "\n".join(lines)


DEFAULT_OUTPUT_PATH = "data/processed/enrichment_results.csv"


if __name__ == "__main__":
    import os

    net = build_network()
    combined = run_all_enrichment_tests(net)
    print(summarize(combined))
    os.makedirs(os.path.dirname(DEFAULT_OUTPUT_PATH), exist_ok=True)
    combined.to_csv(DEFAULT_OUTPUT_PATH, index=False)
    print(f"\nFull results written to {DEFAULT_OUTPUT_PATH}")
