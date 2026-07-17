"""Loads ``data/raw/MuThetaDataSim_m4.txt``, a user-supplied simulated
dataset (source and generating model not yet documented by the project --
see caveats below).

**File-format quirks handled here (discovered by direct inspection, not
assumed):**

1. The file uses bare ``\\r`` (CR-only, classic-Mac-style) line endings, not
   ``\\n`` or ``\\r\\n``. Reading it with a normal line iterator collapses
   the whole file into one line; it must be split on ``\\r`` explicitly.
2. Each record is a run of 16 ``index\\tvalue`` pairs, indices ``0``-``15``
   in order, with no record separator other than the index resetting to
   ``0``. There is no header row and no gene/sample identifier column.
3. Splitting on that reset yields exactly 2418 records for the current
   file (``_m4`` in the filename -- presumably one of several simulation
   scenarios/configurations; siblings such as ``_m1``-``_m3`` may exist but
   are not present in this repo).

**What the two blocks are, and what is/isn't confirmed:**

- Indices 0-9 ("Mu"): 10 continuous values per record, roughly in
  ``[0, 100]``, mean near 50 for most columns. Read as 10 simulated
  feature/measurement values; their real-world meaning is not yet known.
- Indices 10-15 ("Theta"): exactly 6 non-negative values per record that
  sum to ``1.0`` in every record (verified, not assumed) -- unambiguously a
  probability / mixture-weight vector over 6 categories.

**Structural resemblance to this project's real data (suggestive, not
proven):** this project's Connection_Matrix.xlsx defines exactly 6
post-transcriptional "RNA operon" regulators (``REGULATOR_ORDER[-6:]`` in
``network_builder.py``: PostReg5 NCU08295/lhp-1, PostReg6 NCU04504/rrp-3,
PostReg7 NCU03363, PostReg8 NCU04799/pab-1, PostReg9 NCU00919/rok-1,
PostReg10 NCU09349/has-1) -- the same count as the Theta block here, and the
paper's whole method (VTENS) is about estimating soft/ensemble regulator
assignment, which is exactly what a 6-way probability vector would encode.
``compare_theta_dominance_to_real_network()`` below checks this
quantitatively: category 0 of Theta *is* the single dominant category far
more often than a uniform 1/6 (matching lhp-1 being the real network's
dominant post-transcriptional regulator by a wide margin), but the relative
ranking of the other 5 categories does **not** cleanly match the other 5
regulators' real target-count shares. So: the "one dominant category out of
6" structural fact matches; a full one-to-one identity mapping of
theta columns 1-5 to specific regulators is NOT established by this data
alone and should not be assumed.

There is no gene/sample identifier in the raw file, so individual records
cannot yet be joined to ``gene_families.csv``, ``enrichment_results.csv``,
or ``Connection_Matrix.xlsx`` by gene. Row order is preserved on load
(``sim_index`` 0..2417) so that a future gene-name list (same row order)
can be joined on later.
"""
from __future__ import annotations

import os

import pandas as pd

from network_builder import REGULATOR_ORDER

DEFAULT_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "MuThetaDataSim_m4.txt")

N_MU = 10
N_THETA = 6
RECORD_LEN = N_MU + N_THETA  # 16

# The 6 post-transcriptional regulators, in this project's canonical order.
# NOT asserted to be the identity of theta_0..theta_5 -- see module
# docstring caveat. Kept here only so callers can opt into that labeling
# explicitly via `label_theta_as_post_transcriptional=True`.
POST_TRANSCRIPTIONAL_LABELS = REGULATOR_ORDER[-6:]


def _read_records(path: str) -> list[dict[int, float]]:
    with open(path, "r", newline="") as f:
        raw = f.read()

    records: list[dict[int, float]] = []
    current: dict[int, float] = {}
    for line in raw.split("\r"):
        line = line.strip()
        if not line:
            continue
        idx_str, val_str = line.split("\t")
        idx = int(idx_str)
        if idx == 0 and current:
            records.append(current)
            current = {}
        current[idx] = float(val_str)
    if current:
        records.append(current)
    return records


def load_mutheta_sim(
    path: str = DEFAULT_PATH, label_theta_as_post_transcriptional: bool = False
) -> pd.DataFrame:
    """Returns a tidy (n_records x 17) DataFrame: ``mu_0``..``mu_9``,
    ``theta_0``..``theta_5``, and ``dominant_theta`` (argmax column name of
    the theta block for that record).

    If ``label_theta_as_post_transcriptional`` is True, the theta columns
    are named after ``POST_TRANSCRIPTIONAL_LABELS`` instead of
    ``theta_0``..``theta_5`` -- an assumption, not a confirmed mapping (see
    module docstring). Off by default.
    """
    records = _read_records(path)

    rows = []
    for rec in records:
        if sorted(rec) != list(range(RECORD_LEN)):
            raise ValueError(f"Malformed record, expected indices 0-{RECORD_LEN - 1}: {sorted(rec)}")
        rows.append([rec[i] for i in range(RECORD_LEN)])

    mu_cols = [f"mu_{i}" for i in range(N_MU)]
    theta_cols = [f"theta_{i}" for i in range(N_THETA)]
    df = pd.DataFrame(rows, columns=mu_cols + theta_cols)
    df.index.name = "sim_index"

    theta_sums = df[theta_cols].sum(axis=1)
    if not theta_sums.round(6).eq(1.0).all():
        bad = theta_sums[~theta_sums.round(6).eq(1.0)]
        raise ValueError(f"Theta block does not sum to 1.0 for {len(bad)} record(s): {bad.head()}")

    df["dominant_theta"] = df[theta_cols].idxmax(axis=1)

    if label_theta_as_post_transcriptional:
        rename = dict(zip(theta_cols, POST_TRANSCRIPTIONAL_LABELS))
        df = df.rename(columns=rename)
        df["dominant_theta"] = df["dominant_theta"].map(
            dict(zip(theta_cols, POST_TRANSCRIPTIONAL_LABELS))
        )

    return df


def compare_theta_dominance_to_real_network(
    connection_matrix_path: str = "data/raw/Connection_Matrix.xlsx",
) -> pd.DataFrame:
    """Quantitative check of the "theta ~ post-transcriptional regulator
    assignment" hypothesis: how often each theta column wins the argmax in
    the simulated data, vs. each real post-transcriptional regulator's
    share of real target genes (from Connection_Matrix.xlsx). Positional
    only (theta_0 vs POST_TRANSCRIPTIONAL_LABELS[0], etc.) -- see module
    docstring for why this pairing is not confirmed.
    """
    from data_loader import REGULATOR_INFO, load_clean

    sim = load_mutheta_sim()
    theta_cols = [f"theta_{i}" for i in range(N_THETA)]
    sim_share = sim["dominant_theta"].value_counts(normalize=True).reindex(theta_cols).fillna(0.0)

    targets, _, _ = load_clean(connection_matrix_path)
    targets = targets.reindex(REGULATOR_ORDER)
    post_counts = targets.loc[POST_TRANSCRIPTIONAL_LABELS].sum(axis=1)
    real_share = post_counts / post_counts.sum()

    return pd.DataFrame(
        {
            "theta_column": theta_cols,
            "candidate_regulator": POST_TRANSCRIPTIONAL_LABELS,
            "sim_dominant_share": sim_share.values,
            "real_target_share": real_share.values,
        }
    ).set_index("theta_column")


if __name__ == "__main__":
    df = load_mutheta_sim()
    print(f"Loaded {len(df)} records")
    print(df.head())
    print()
    print("Dominant theta column distribution:")
    print(df["dominant_theta"].value_counts(normalize=True).round(3))
    print()
    print("Comparison to real post-transcriptional regulator target shares:")
    print(compare_theta_dominance_to_real_network().round(3))
