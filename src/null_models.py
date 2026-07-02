"""Degree-preserving null model for the observed FFL count.

Novel extension beyond the paper (see README): the paper reports 71 FFLs
(here: an observed count from `ffl_analysis.count_ffls`) but never asks
whether that count is more than expected by chance given each regulator's
target count. This module answers that with a configuration-model-style
null: repeatedly apply degree-preserving double-edge-swaps to T (holding R,
the fixed WCC -> regulator structure, untouched) and recompute the FFL
count on each randomized network to build an empirical null distribution.

Double-edge-swap on a binary regulator x target matrix, applied to a pair
of rows (r1, r2) and a pair of columns (c1, c2) with T[r1,c1]=1, T[r2,c1]=0,
T[r2,c2]=1, T[r1,c2]=0 -> flip all four cells. This exactly preserves every
row sum (each regulator's target count) and every column sum (each gene's
in-degree, i.e. how many regulators target it), which is the standard
degree-preserving randomization for bipartite/binary adjacency matrices
(also known as the "curveball"/double-edge-swap algorithm).
"""
from __future__ import annotations

import os
import time
from dataclasses import dataclass

import numpy as np

from network_builder import ClockNetwork, WCC_INDEX, build_network

DEFAULT_CACHE_PATH = "data/processed/null_model_results.npz"


@dataclass
class NullModelResult:
    observed_ffl_count: int
    null_counts: np.ndarray       # (n_iterations,) FFL count per randomized network
    n_iterations: int
    n_swaps_per_iteration: int
    mean_null: float
    std_null: float
    z_score: float
    p_value_upper: float          # P(null >= observed), empirical, one-sided
    p_value_two_sided: float      # 2 * min(P(null>=observed), P(null<=observed)), capped at 1
    elapsed_seconds: float


def _double_edge_swap(T: np.ndarray, n_swaps: int, rng: np.random.Generator) -> int:
    """In-place degree-preserving double-edge-swap on boolean matrix T.

    Returns the number of successful swaps (attempts can fail when the
    randomly chosen row pair has no valid column pair to swap, e.g. one
    row's target set is a subset of the other's).
    """
    n_rows = T.shape[0]
    successes = 0
    attempts = 0
    max_attempts = n_swaps * 20  # generous ceiling so we don't loop forever on tiny/degenerate inputs
    while successes < n_swaps and attempts < max_attempts:
        attempts += 1
        r1, r2 = rng.choice(n_rows, size=2, replace=False)
        row1, row2 = T[r1], T[r2]
        diff1 = np.flatnonzero(row1 & ~row2)
        if diff1.size == 0:
            continue
        diff2 = np.flatnonzero(row2 & ~row1)
        if diff2.size == 0:
            continue
        c1 = diff1[rng.integers(diff1.size)]
        c2 = diff2[rng.integers(diff2.size)]
        T[r1, c1], T[r2, c1] = False, True
        T[r2, c2], T[r1, c2] = False, True
        successes += 1
    return successes


def _ffl_count_from_matrix(T: np.ndarray, wcc_index: int) -> int:
    """Fast FFL count: sum over non-WCC rows r of |targets(WCC) & targets(r)|.

    Equivalent to summing the Hadamard product (R @ T) * T used in
    ffl_analysis.count_ffls, since R's only nonzero row is WCC's (see that
    module's docstring for the derivation).
    """
    wcc_row = T[wcc_index]
    total = 0
    for r in range(T.shape[0]):
        if r == wcc_index:
            continue
        total += int(np.count_nonzero(wcc_row & T[r]))
    return total


def run_null_model(
    net: ClockNetwork,
    n_iterations: int = 1000,
    n_swaps_per_iteration: int | None = None,
    swap_multiplier: float = 5.0,
    seed: int | None = 0,
    verbose: bool = False,
) -> NullModelResult:
    """Build an empirical null distribution of FFL counts via degree-preserving rewiring.

    Parameters
    ----------
    n_iterations : number of independent randomized networks to generate.
    n_swaps_per_iteration : explicit swap count per randomized network. If
        None, defaults to ``swap_multiplier * total_edges`` (a standard
        rule of thumb for adequate mixing in double-edge-swap null models).
    swap_multiplier : used only when n_swaps_per_iteration is None.
    seed : seed for the numpy Generator (reproducibility).
    """
    T0 = net.T.values.astype(bool)
    n_edges = int(T0.sum())
    if n_swaps_per_iteration is None:
        n_swaps_per_iteration = int(round(swap_multiplier * n_edges))

    observed = _ffl_count_from_matrix(T0, WCC_INDEX)

    rng = np.random.default_rng(seed)
    null_counts = np.empty(n_iterations, dtype=int)

    start = time.time()
    for i in range(n_iterations):
        T = T0.copy()
        _double_edge_swap(T, n_swaps_per_iteration, rng)
        null_counts[i] = _ffl_count_from_matrix(T, WCC_INDEX)
        if verbose and (i + 1) % max(1, n_iterations // 10) == 0:
            elapsed = time.time() - start
            print(f"  [{i + 1}/{n_iterations}] elapsed={elapsed:.1f}s")
    elapsed = time.time() - start

    mean_null = float(null_counts.mean())
    std_null = float(null_counts.std(ddof=1))
    z = (observed - mean_null) / std_null if std_null > 0 else float("nan")

    p_upper = float((null_counts >= observed).sum() + 1) / (n_iterations + 1)
    p_lower = float((null_counts <= observed).sum() + 1) / (n_iterations + 1)
    p_two_sided = min(1.0, 2 * min(p_upper, p_lower))

    return NullModelResult(
        observed_ffl_count=observed,
        null_counts=null_counts,
        n_iterations=n_iterations,
        n_swaps_per_iteration=n_swaps_per_iteration,
        mean_null=mean_null,
        std_null=std_null,
        z_score=z,
        p_value_upper=p_upper,
        p_value_two_sided=p_two_sided,
        elapsed_seconds=elapsed,
    )


def save_result(result: NullModelResult, path: str = DEFAULT_CACHE_PATH) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(
        path,
        observed_ffl_count=result.observed_ffl_count,
        null_counts=result.null_counts,
        n_iterations=result.n_iterations,
        n_swaps_per_iteration=result.n_swaps_per_iteration,
        mean_null=result.mean_null,
        std_null=result.std_null,
        z_score=result.z_score,
        p_value_upper=result.p_value_upper,
        p_value_two_sided=result.p_value_two_sided,
        elapsed_seconds=result.elapsed_seconds,
    )


def load_result(path: str = DEFAULT_CACHE_PATH) -> NullModelResult | None:
    if not os.path.exists(path):
        return None
    data = np.load(path)
    return NullModelResult(
        observed_ffl_count=int(data["observed_ffl_count"]),
        null_counts=data["null_counts"],
        n_iterations=int(data["n_iterations"]),
        n_swaps_per_iteration=int(data["n_swaps_per_iteration"]),
        mean_null=float(data["mean_null"]),
        std_null=float(data["std_null"]),
        z_score=float(data["z_score"]),
        p_value_upper=float(data["p_value_upper"]),
        p_value_two_sided=float(data["p_value_two_sided"]),
        elapsed_seconds=float(data["elapsed_seconds"]),
    )


def summarize(result: NullModelResult) -> str:
    lines = [
        f"Degree-preserving null model ({result.n_iterations} randomized networks, "
        f"{result.n_swaps_per_iteration} swaps each, {result.elapsed_seconds:.1f}s)",
        f"Observed FFL count: {result.observed_ffl_count}",
        f"Null distribution: mean={result.mean_null:.2f}, std={result.std_null:.2f}",
        f"Z-score: {result.z_score:.3f}",
        f"Empirical p-value (upper-tail, P(null >= observed)): {result.p_value_upper:.4g}",
        f"Empirical p-value (two-sided): {result.p_value_two_sided:.4g}",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    net = build_network()
    result = run_null_model(net, n_iterations=1000, verbose=True)
    print()
    print(summarize(result))
    save_result(result)
    print(f"\nSaved to {DEFAULT_CACHE_PATH}")
