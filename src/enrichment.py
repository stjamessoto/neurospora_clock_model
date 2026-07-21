"""Pairwise, triple-wise, and quartet-wise regulator co-targeting enrichment tests.

Novel extension beyond the paper (see README): the paper never tested
whether any two, three, or four regulators share more target genes than
expected by chance. Null model: each regulator's target set is treated as
a fixed-size random subset of the N-gene universe, drawn independently of
every other regulator's target set (their sizes -- "row sums" of T -- are
taken as given, exactly as in the paper's fitted network).

- Pairwise (55 = C(11,2) tests): exact hypergeometric test (equivalent to
  a two-sided Fisher's exact test on the 2x2 contingency table for each
  regulator pair), via scipy.stats.hypergeom.
- Triple-wise (165 = C(11,3) tests) and quartet-wise (330 = C(11,4) tests):
  exact distribution of a k-way intersection size under the same
  independent-random-subset null, computed by sequentially convolving
  hypergeometric laws:
      P(|A1 n A2| = m)                 ~ Hypergeom(N, n1, n2)
      P(|A3 n (A1 n A2)| = k | m)      ~ Hypergeom(N, m, n3)
      P(|A1 n A2 n A3| = k) = sum_m P(|A1 n A2|=m) * P(k | m)
  and, for a fourth set, repeating the same conditioning step against the
  running (A1 n A2 n A3) distribution. This is exact (not a normal/
  simulation approximation) because the sets are independent draws, so the
  "population" being intersected with each next set is exactly the
  (random-sized) running intersection.

All 55 + 165 + 330 = 550 p-values (two-sided) are pooled and Benjamini-
Hochberg FDR corrected at alpha = 0.05 (statsmodels multipletests),
mirroring the significance-table format used in the companion
S. cerevisiae subnetwork-motif project.
"""
from __future__ import annotations

from itertools import combinations
from typing import Sequence

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


def _nway_intersection_pmf(N: int, sizes: Sequence[int]) -> np.ndarray:
    """Exact pmf of |A_1 n A_2 n ... n A_k| for independent random subsets of
    the given sizes drawn from a population of N.

    Starts from the exact pairwise distribution of |A_1 n A_2|, then folds
    in each remaining set by conditioning on the running intersection size
    (see module docstring). Returns `pmf` where pmf[k] = P(|intersection| = k).
    """
    n1, n2, *rest = sizes
    m_min = max(0, n1 + n2 - N)
    m_max = min(n1, n2)
    pmf = np.zeros(m_max + 1)
    pmf[m_min:] = hypergeom.pmf(np.arange(m_min, m_max + 1), N, n1, n2)

    for n_next in rest:
        max_k = min(len(pmf) - 1, n_next)
        new_pmf = np.zeros(max_k + 1)
        for m, pm in enumerate(pmf):
            if pm <= 0:
                continue
            k_max_m = min(m, n_next)
            ks = np.arange(0, k_max_m + 1)
            p_k_given_m = hypergeom.pmf(ks, N, m, n_next)
            new_pmf[: k_max_m + 1] += pm * p_k_given_m
        pmf = new_pmf
    return pmf


def _nway_pvalue(observed: int, N: int, sizes: Sequence[int]) -> tuple[float, float, float]:
    """Two-sided exact p-value for a k-way overlap (k >= 3), plus its null expectation.

    Returns (expected, p_upper, p_two_sided).
    """
    expected = float(np.prod(sizes)) / (N ** (len(sizes) - 1))
    pmf = _nway_intersection_pmf(N, sizes)
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
                "reg_d": None,
                "n_a": n1,
                "n_b": n2,
                "n_c": np.nan,
                "n_d": np.nan,
                "observed_overlap": observed,
                "expected_overlap": expected,
                "direction": "enriched" if observed >= expected else "depleted",
                "p_upper": p_upper,
                "p_value": p_two_sided,
            }
        )
    return pd.DataFrame(rows)


def _kway_enrichment(net: ClockNetwork, k: int, test_type: str) -> pd.DataFrame:
    """Shared implementation for triple (k=3) and quartet (k=4) enrichment."""
    T = net.T
    N = T.shape[1]
    labels = list(T.index)
    names = {lab: REGULATOR_INFO[lab]["gene"] for lab in labels}
    sizes = T.sum(axis=1)
    bool_sets = {lab: set(T.columns[T.loc[lab].values.astype(bool)]) for lab in labels}

    rows = []
    for combo in combinations(labels, k):
        observed = len(set.intersection(*(bool_sets[lab] for lab in combo)))
        combo_sizes = [int(sizes[lab]) for lab in combo]
        expected, p_upper, p_two_sided = _nway_pvalue(observed, N, combo_sizes)
        combo_padded = list(combo) + [None] * (4 - k)
        sizes_padded = combo_sizes + [np.nan] * (4 - k)
        rows.append(
            {
                "test_type": test_type,
                "regulators": " + ".join(names[lab] for lab in combo),
                "reg_a": combo_padded[0],
                "reg_b": combo_padded[1],
                "reg_c": combo_padded[2],
                "reg_d": combo_padded[3],
                "n_a": sizes_padded[0],
                "n_b": sizes_padded[1],
                "n_c": sizes_padded[2],
                "n_d": sizes_padded[3],
                "observed_overlap": observed,
                "expected_overlap": expected,
                "direction": "enriched" if observed >= expected else "depleted",
                "p_upper": p_upper,
                "p_value": p_two_sided,
            }
        )
    return pd.DataFrame(rows)


def triple_enrichment(net: ClockNetwork) -> pd.DataFrame:
    return _kway_enrichment(net, 3, "triple")


def quartet_enrichment(net: ClockNetwork) -> pd.DataFrame:
    return _kway_enrichment(net, 4, "quartet")


def run_all_enrichment_tests(net: ClockNetwork, alpha: float = ALPHA) -> pd.DataFrame:
    """Run all 55 pairwise + 165 triple + 330 quartet tests, pool, and apply BH-FDR correction."""
    pairs = pairwise_enrichment(net)
    triples = triple_enrichment(net)
    quartets = quartet_enrichment(net)
    combined = pd.concat([pairs, triples, quartets], ignore_index=True)

    assert len(combined) == 55 + 165 + 330, f"expected 550 tests, got {len(combined)}"

    reject, q_values, _, _ = multipletests(combined["p_value"].values, alpha=alpha, method="fdr_bh")
    combined["q_value"] = q_values
    combined["significant_bh"] = reject
    combined = combined.sort_values("q_value", ascending=True).reset_index(drop=True)
    return combined


def summarize(combined: pd.DataFrame, alpha: float = ALPHA) -> str:
    n_sig = int(combined["significant_bh"].sum())
    counts = combined["test_type"].value_counts()
    sig_counts = combined.loc[combined["significant_bh"], "test_type"].value_counts()

    lines = [
        f"Ran {len(combined)} tests "
        f"({counts.get('pair', 0)} pairs + {counts.get('triple', 0)} triples + "
        f"{counts.get('quartet', 0)} quartets), BH-FDR corrected at alpha={alpha}.",
        f"Significant after correction: {n_sig} total "
        f"({sig_counts.get('pair', 0)} pairs, {sig_counts.get('triple', 0)} triples, "
        f"{sig_counts.get('quartet', 0)} quartets).",
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
