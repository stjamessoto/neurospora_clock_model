"""Loads the binding-strength-by-function table
(``data/raw/binding_strength_by_function.csv``, user-supplied, transcribed
from the paper's supplementary "Binding strength table") -- 21 KEGG/GO
functional categories x 11 regulators, continuous values in ``[0, 1)``.

**What these numbers are, and are not.** Al-Omari et al.'s Variable Topology
Ensemble Method (VTENS) does not fit each regulator-target binding strength
(mu_k) as an independent parameter with its own marginal posterior later
tested against zero. Under VTENS's "supernet" formulation, every regulator
is initially allowed to regulate every candidate target/function, and the
network topology and kinetic parameters are estimated *jointly* by MCMC --
each accepted sample is one whole candidate topology consistent with the
mRNA time-course data. The binding strengths reported here are the fitted
ensemble's joint point estimates (used afterward to pick out the dominant
regulator for each gene/function and build `Connection_Matrix.xlsx`), not
independent per-parameter distributions -- there is no marginal ensemble
distribution per (regulator, function) cell to test for excluding zero, and
this table should never be read as the result of such a test. A `0.0000` in
this table means "not the estimated dominant regulator for that function,"
not "ensemble excluded zero at some confidence level."
"""
from __future__ import annotations

import os

import pandas as pd

from data_loader import REGULATOR_INFO
from network_builder import REGULATOR_ORDER

DEFAULT_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "binding_strength_by_function.csv")

# The raw file's columns are NCU IDs / "WCC"; map them to this project's
# canonical regulator labels (REGULATOR_ORDER), the same mapping
# data_loader.py uses for the connection matrix.
_NCU_TO_LABEL = {info["ncu"]: label for label, info in REGULATOR_INFO.items() if info["ncu"]}
_NCU_TO_LABEL["WCC"] = "Wcc"


def load_binding_strength_by_function(path: str = DEFAULT_PATH) -> pd.DataFrame:
    """Returns a (n_functions x 11 regulators) DataFrame of continuous
    binding-strength point estimates -- columns renamed to canonical
    regulator labels and reordered to ``REGULATOR_ORDER``, index is the
    function/pathway name."""
    df = pd.read_csv(path, index_col="Function")
    df = df.rename(columns=_NCU_TO_LABEL)
    return df[REGULATOR_ORDER]


if __name__ == "__main__":
    table = load_binding_strength_by_function()
    print(f"Loaded {table.shape[0]} functions x {table.shape[1]} regulators")
    print(table)
