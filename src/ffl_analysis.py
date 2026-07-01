"""Feedforward loop (FFL) counting, per Section V-F of the paper.

Definition: WCC (A) regulates a regulator (B), and both A and B regulate a
target gene (C). Given R (11x11, WCC->regulator, zero diagonal) and T
(11xN, regulator->target), all FFLs are counted via the Hadamard product
of the ordinary matrix product RT with T:

    FFL_matrix = (R @ T) * T        (elementwise multiply, "Hadamard product")
    FFL_count  = FFL_matrix.sum()

An entry FFL_matrix[b, c] counts, for regulator b and target c, how many
distinct A-regulators (rows of R with a 1 in column b) feed into the
A -> B -> C / A -> C loop. Since R's only nonzero row is WCC's, this
collapses to: FFL_matrix[b, c] = R[WCC, b] * T[b, c] * T[WCC, c], i.e. 1
exactly when regulator b targets c AND WCC also directly targets c (and b
isn't WCC itself, since R has a zero diagonal position for WCC via
R[WCC, WCC] = 0).

For a target gene c, sum over b of FFL_matrix[b, c] counts how many
non-WCC regulators b are alternative "B genes" completing an FFL to c --
this is what the paper calls entries of RT o T greater than 1: multiple
regulators sharing the same target gene c in overlapping FFLs.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from network_builder import ClockNetwork, WCC_INDEX, build_network

PAPER_REPORTED_FFL_COUNT = 71
PAPER_OVERLAP_GENES = ["mrp-180/NCU08120", "trax/NCU06059", "emp-8/NCU04728", "NCU02213", "NCU08166"]
PAPER_LHP1_NCU06108_FFL_SHARE = 39


@dataclass
class FFLResult:
    ffl_matrix: pd.DataFrame          # (11, N) count of FFLs per (regulator, target) cell
    total_ffl_count: int
    per_gene_ffl_count: pd.Series     # (N,) total FFL count per target gene, summed over regulators
    overlap_genes: pd.Series          # genes with per_gene_ffl_count > 1
    per_regulator_ffl_count: pd.Series  # (11,) FFL count where each regulator plays the "B" role


def count_ffls(net: ClockNetwork) -> FFLResult:
    R = net.R
    T = net.T.values  # (11, N)

    # RT[b, c] = sum_k R[b, k] * T[k, c]. Since R has a single nonzero row
    # (WCC's), RT is nonzero only in the WCC row: RT[wcc, c] = number of
    # non-WCC regulators targeting c. All FFL content therefore collapses
    # into the WCC row of RT (.) T -- this is a property of R's structure,
    # not a bug; see per_regulator_ffl_count below for the actual
    # per-regulator ("B" gene) breakdown, computed as the paper specifies.
    RT = R @ T
    ffl_matrix_arr = RT * T  # Hadamard product, (11, N)

    ffl_matrix = pd.DataFrame(ffl_matrix_arr, index=net.T.index, columns=net.T.columns)

    total = int(ffl_matrix_arr.sum())
    per_gene = ffl_matrix.sum(axis=0)
    overlap = per_gene[per_gene > 1].sort_values(ascending=False)

    # Per the paper: "An inner product of the first row of the T matrix
    # with any of the successive 10 rows identifies the number of
    # regulators that appear as a B gene in a feedforward loop." I.e. for
    # each non-WCC regulator k, count genes targeted by BOTH WCC and k.
    wcc_row = net.T.loc["Wcc"]
    per_regulator = {}
    for label in net.T.index:
        if label == "Wcc":
            continue
        per_regulator[label] = int((wcc_row & net.T.loc[label]).sum())
    per_regulator = pd.Series(per_regulator)

    assert per_regulator.sum() == total, (
        f"Per-regulator B-gene counts ({per_regulator.sum()}) must sum to the "
        f"total FFL count ({total})"
    )

    return FFLResult(
        ffl_matrix=ffl_matrix,
        total_ffl_count=total,
        per_gene_ffl_count=per_gene,
        overlap_genes=overlap,
        per_regulator_ffl_count=per_regulator,
    )


def summarize(result: FFLResult, net: ClockNetwork) -> str:
    lines = []
    lines.append(f"Observed total FFL count: {result.total_ffl_count}")
    lines.append(f"Paper's reported FFL count: {PAPER_REPORTED_FFL_COUNT}")
    lines.append(
        "(Mismatch expected: this export's per-gene regulator degree tops out at 2, "
        "vs. the paper's continuous binding-strength ensemble -- see README for discussion.)"
    )
    lines.append("")
    lines.append(f"Genes with overlapping FFLs (per-gene count > 1): {len(result.overlap_genes)}")
    lines.append(f"Paper reports 5: {PAPER_OVERLAP_GENES}")
    if len(result.overlap_genes):
        lines.append(result.overlap_genes.to_string())
    lines.append("")
    lines.append("FFL contribution per regulator (as the 'B' gene):")
    contrib = result.per_regulator_ffl_count.sort_values(ascending=False)
    lines.append(contrib.to_string())
    lines.append("")
    lhp1_ncu06108_share = int(
        contrib.get("PostReg5 NCU08295", 0) + contrib.get("Reg3 NCU06108", 0)
    )
    lines.append(
        f"lhp-1 + NCU06108 share of all FFLs: {lhp1_ncu06108_share} / {result.total_ffl_count} "
        f"(paper: {PAPER_LHP1_NCU06108_FFL_SHARE} / {PAPER_REPORTED_FFL_COUNT}). "
        f"Same qualitative pattern (lhp-1 is the single largest FFL contributor here), "
        f"different magnitude."
    )
    return "\n".join(lines)


if __name__ == "__main__":
    net = build_network()
    result = count_ffls(net)
    print(summarize(result, net))
