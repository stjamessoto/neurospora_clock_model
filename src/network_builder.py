"""Construct the R (WCC -> regulator) and T (regulator -> target) matrices.

Per the paper (Section III, Section V-F): all 11 non-WCC regulators are
activated by the master regulator WCC with fixed binding constant mu_0 = 1,
and "the regulatory CCG proteins do not regulate each other." So R is
trivial and is NOT present as data in Connection_Matrix.xlsx -- it is
constructed here directly from that stated topology, not inferred.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from data_loader import REGULATOR_INFO, load_clean

REGULATOR_ORDER = list(REGULATOR_INFO.keys())  # Wcc first, then Reg1..Reg4, PostReg5..PostReg10
WCC_INDEX = REGULATOR_ORDER.index("Wcc")


@dataclass
class ClockNetwork:
    R: np.ndarray          # (11, 11) WCC -> regulator adjacency, zero diagonal
    T: pd.DataFrame        # (11, N) regulator -> target adjacency, binary
    core_clock: pd.DataFrame  # (11, 3) regulator -> core clock gene adjacency
    regulator_labels: list[str]   # raw row labels, e.g. "PostReg5 NCU08295"
    regulator_names: list[str]    # gene names, e.g. "lhp-1"
    regulator_class: list[str]    # "master" / "transcriptional" / "post-transcriptional"


def build_R(n_regulators: int = 11, wcc_index: int = WCC_INDEX) -> np.ndarray:
    """WCC (row/col `wcc_index`) activates every other regulator; zero diagonal.

    Regulators do not regulate each other, per the paper.
    """
    R = np.zeros((n_regulators, n_regulators), dtype=int)
    R[wcc_index, :] = 1
    R[:, wcc_index] = 0  # WCC is not self-activated or activated by others
    np.fill_diagonal(R, 0)
    return R


def build_network(path: str = "data/raw/Connection_Matrix.xlsx") -> ClockNetwork:
    targets, core_clock, report = load_clean(path)

    # Reorder rows to a canonical order (Wcc, Reg1..Reg4, PostReg5..PostReg10)
    targets = targets.reindex(REGULATOR_ORDER)
    core_clock = core_clock.reindex(REGULATOR_ORDER)

    R = build_R(n_regulators=len(REGULATOR_ORDER), wcc_index=WCC_INDEX)

    names = [REGULATOR_INFO[r]["gene"] for r in REGULATOR_ORDER]
    classes = [REGULATOR_INFO[r]["class"] for r in REGULATOR_ORDER]

    return ClockNetwork(
        R=R,
        T=targets,
        core_clock=core_clock,
        regulator_labels=REGULATOR_ORDER,
        regulator_names=names,
        regulator_class=classes,
    )


if __name__ == "__main__":
    net = build_network()
    print("R matrix (WCC -> regulators):")
    print(pd.DataFrame(net.R, index=net.regulator_names, columns=net.regulator_names))
    print()
    print("T matrix shape:", net.T.shape)
    print("Regulator target counts:")
    print(pd.Series(net.T.sum(axis=1).values, index=net.regulator_names))
