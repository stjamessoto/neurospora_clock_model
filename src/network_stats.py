"""Descriptive network statistics: degree distribution, hub sizes, and the
singly- vs multiply-regulated gene split. Compare to Section V-E ("The View
of the Clock Network in Toto With Its 12 Hubs") of the paper.

Paper's reported numbers for this section:
- Most connected node: lhp-1 (NCU08295), 768 connections.
- Least connected regulatory node: rok-1 (NCU00919), 54 connections.
- Average connections per node: 2.23.
- Degree distribution: long right tail, power-law exponent approx -1.4.
- 2686 genes regulated by exactly one regulator; 694 by >=2 (2686+694=3380).

As documented in the README, this export's per-gene regulator degree tops
out at 2 (vs. the paper's up-to-11), so the degree distribution here is
structurally much less rich than the paper's -- comparisons below are
reported, not forced to match.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from network_builder import ClockNetwork, build_network

PAPER_MAX_HUB = 768          # lhp-1
PAPER_MIN_HUB = 54           # rok-1
PAPER_AVG_NODE_DEGREE = 2.23
PAPER_POWER_LAW_EXPONENT = -1.4
PAPER_SINGLY_REGULATED = 2686
PAPER_MULTIPLY_REGULATED = 694
PAPER_TOTAL_GENES = 3380


@dataclass
class NetworkStatsResult:
    hub_sizes: pd.Series             # regulator gene name -> target count (row sums), sorted desc
    gene_degrees: pd.Series          # target gene id -> number of regulators targeting it (col sums)
    singly_regulated: int
    multiply_regulated: int
    degree_distribution: pd.Series   # degree value -> count of nodes (genes + regulators) at that degree
    power_law_exponent: float
    power_law_intercept: float
    power_law_r_squared: float
    average_node_degree: float


def compute_network_stats(net: ClockNetwork) -> NetworkStatsResult:
    T = net.T
    hub_sizes = T.sum(axis=1)
    hub_sizes.index = net.regulator_names
    hub_sizes = hub_sizes.sort_values(ascending=False)

    gene_degrees = T.sum(axis=0)

    singly = int((gene_degrees == 1).sum())
    multiply = int((gene_degrees >= 2).sum())

    # Combined degree distribution over all nodes in the bipartite graph:
    # 11 regulator nodes (degree = target count) + N target-gene nodes
    # (degree = number of regulators targeting them), as in the paper's
    # Figure 11 degree distribution.
    all_degrees = pd.concat(
        [hub_sizes.reset_index(drop=True), gene_degrees.reset_index(drop=True)],
        ignore_index=True,
    )
    degree_counts = all_degrees.value_counts().sort_index()

    mask = degree_counts.index > 0
    log_k = np.log10(degree_counts.index[mask].to_numpy(dtype=float))
    log_freq = np.log10(degree_counts.values[mask].astype(float))
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_k, log_freq)

    avg_degree = float(all_degrees.mean())

    return NetworkStatsResult(
        hub_sizes=hub_sizes,
        gene_degrees=gene_degrees,
        singly_regulated=singly,
        multiply_regulated=multiply,
        degree_distribution=degree_counts,
        power_law_exponent=float(slope),
        power_law_intercept=float(intercept),
        power_law_r_squared=float(r_value ** 2),
        average_node_degree=avg_degree,
    )


def summarize(result: NetworkStatsResult) -> str:
    lines = [
        "Hub sizes (regulator -> target count), descending:",
        result.hub_sizes.to_string(),
        "",
        f"Largest hub here: {result.hub_sizes.index[0]} ({int(result.hub_sizes.iloc[0])}) "
        f"-- paper: lhp-1 ({PAPER_MAX_HUB})",
        f"Smallest hub here: {result.hub_sizes.index[-1]} ({int(result.hub_sizes.iloc[-1])}) "
        f"-- paper: rok-1 ({PAPER_MIN_HUB})",
        "",
        f"Singly-regulated genes: {result.singly_regulated} -- paper: {PAPER_SINGLY_REGULATED}",
        f"Multiply-regulated genes (>=2 regulators): {result.multiply_regulated} "
        f"-- paper: {PAPER_MULTIPLY_REGULATED}",
        f"Total genes: {result.singly_regulated + result.multiply_regulated} "
        f"-- paper: {PAPER_TOTAL_GENES}",
        "",
        f"Average node degree (regulators + genes): {result.average_node_degree:.3f} "
        f"-- paper: {PAPER_AVG_NODE_DEGREE}",
        f"Fitted power-law exponent (log-log linreg over degree value frequencies): "
        f"{result.power_law_exponent:.3f} (R^2={result.power_law_r_squared:.3f}) "
        f"-- paper: {PAPER_POWER_LAW_EXPONENT}",
        "",
        "Caveat: this export's per-gene regulator degree only takes values {1, 2} "
        "(see README), so the fitted exponent above is estimated from a handful of "
        "distinct degree values (2 gene-side + up to 11 regulator-side) rather than a "
        "rich continuous tail -- treat it as illustrative, not a robust power-law fit.",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    net = build_network()
    result = compute_network_stats(net)
    print(summarize(result))
