"""Loads ``data/raw/MuThetaDataSim_m4.txt`` (a user-supplied simulated
dataset) and ``data/raw/ALLNCU.txt`` (the gene-identifier list for it,
supplied later, same row order), then tests -- not just gestures at -- the
hypothesis that Theta encodes per-gene post-transcriptional regulator
assignment.

**File-format quirks handled here (discovered by direct inspection, not
assumed):**

1. ``MuThetaDataSim_m4.txt`` uses bare ``\\r`` (CR-only, classic-Mac-style)
   line endings, not ``\\n`` or ``\\r\\n``. Reading it with a normal line
   iterator collapses the whole file into one line; it must be split on
   ``\\r`` explicitly.
2. Each record is a run of 16 ``index\\tvalue`` pairs, indices ``0``-``15``
   in order, with no record separator other than the index resetting to
   ``0``. There is no header row and no gene/sample identifier column in
   this file itself.
3. Splitting on that reset yields exactly 2418 records for the current
   file (``_m4`` in the filename -- presumably one of several simulation
   scenarios/configurations; siblings such as ``_m1``-``_m3`` may exist but
   are not present in this repo).
4. ``ALLNCU.txt`` is a plain newline-delimited list of NCU gene IDs, also
   exactly 2418 lines. Row order is a one-to-one, per-gene correspondence
   with ``MuThetaDataSim_m4.txt`` -- confirmed by the data source (not just
   inferred from the matching line counts) to be the indexing used
   throughout the analysis to associate each Mu/Theta record with its
   gene (see ``load_mutheta_sim_with_genes``). It has 202 duplicate lines
   (2216 distinct genes across 2418 rows) -- some genes appear at more than
   one ``sim_index``, each with its own Mu/Theta values.

**What the two blocks are, and what is/isn't confirmed:**

- Indices 0-9 ("Mu"): 10 continuous values per record, roughly in
  ``[0, 100]``. Per the data source: these are **conditional ensemble
  averages**, not unconditional means, over an ensemble of 2000 MCMC runs.
  For each ensemble member, the run's code finds the regulator with the
  largest inferred binding strength for that gene; a given Mu column's
  value for the gene is the average of that regulator's binding strength
  taken *only* over the ensemble members in which it was the dominant
  regulator for that gene, not an average over all 2000 members. This
  confirms what the values represent in general (per-gene, per-regulator
  conditional-dominance binding strengths, joined to genes via the
  ``ALLNCU.txt`` row order above) -- what is still NOT confirmed is which
  of this project's 10 non-WCC regulators each of ``mu_0``..``mu_9``
  positionally corresponds to; treat the column-to-regulator identity as
  open even though the column's *meaning* is now known.
- Indices 10-15 ("Theta"): exactly 6 non-negative values per record that
  sum to ``1.0`` in every record (verified, not assumed) -- unambiguously a
  probability / mixture-weight vector over 6 categories. The data source's
  clarification above was specific to Mu; Theta's semantic identity is
  still tested (and refuted, for the tested hypothesis) below.

**Confirmed by the gene list:** every one of the 2216 distinct genes in
``ALLNCU.txt`` is present as a real target-gene column in the cleaned
``Connection_Matrix.xlsx`` matrix (100% overlap) -- strong evidence this
dataset is drawn from the same gene universe as the rest of this project,
not an unrelated import.

**Tested and NOT confirmed (this matters -- read before using Theta as a
regulator label):** this project's Connection_Matrix.xlsx defines exactly 6
post-transcriptional "RNA operon" regulators (``REGULATOR_ORDER[-6:]`` in
``network_builder.py``: PostReg5 NCU08295/lhp-1, PostReg6 NCU04504/rrp-3,
PostReg7 NCU03363, PostReg8 NCU04799/pab-1, PostReg9 NCU00919/rok-1,
PostReg10 NCU09349/has-1) -- the same count as the Theta block, which
initially looked like a promising match (theta_0 dominates the argmax at a
rate that echoes lhp-1's real dominance). ``best_theta_to_regulator_mapping()``
tests this directly, per gene, against the real data: for the 887 genes in
the sim's gene list that are targeted by exactly one of the 6 real
post-transcriptional regulators, it brute-forces all 720 permutations of
"theta column -> regulator label" and scores each by how often
argmax(mean Theta for that gene) matches that gene's real regulator.

**Result: the best possible permutation scores ~0.29 accuracy -- worse than
the trivial baseline of always guessing lhp-1 (~0.56, since lhp-1 alone
targets 500/887 = 56% of these genes), and only marginally better than the
mean-over-all-permutations of 1/6 (chance).** This refutes the "Theta is a
per-gene post-transcriptional-regulator assignment posterior" hypothesis --
the shared category count and the single-dominant-category shape were
apparently coincidental, not a real semantic match. What Mu and Theta
actually represent is still an open question; per-gene real regulator
identity is *not* it, at least not directly. Flagged rather than quietly
dropped, per this project's convention of reporting discrepancies instead
of forcing numbers to match.
"""
from __future__ import annotations

import itertools
import os

import pandas as pd

from network_builder import REGULATOR_ORDER

DEFAULT_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "MuThetaDataSim_m4.txt")
GENE_LIST_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "ALLNCU.txt")

N_MU = 10
N_THETA = 6
RECORD_LEN = N_MU + N_THETA  # 16

# The 6 post-transcriptional regulators, in this project's canonical order.
# NOT confirmed to be the identity of theta_0..theta_5 -- see module
# docstring; best_theta_to_regulator_mapping() tests and refutes the
# straightforward positional version of this assumption. Kept here so
# callers can still opt into that (unconfirmed) labeling explicitly via
# `label_theta_as_post_transcriptional=True`.
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


def load_gene_list(path: str = GENE_LIST_PATH) -> list[str]:
    """Returns the NCU gene IDs in ``ALLNCU.txt``, in file order (duplicates
    kept as-is -- some genes appear at more than one sim row)."""
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def load_mutheta_sim_with_genes(
    mutheta_path: str = DEFAULT_PATH, gene_list_path: str = GENE_LIST_PATH
) -> pd.DataFrame:
    """``load_mutheta_sim()`` with a ``gene`` column joined on by row order.

    Raises if the two files don't have matching lengths -- the join is
    positional (``sim_index`` i <-> line i of ``ALLNCU.txt``), so a length
    mismatch means the correspondence can't be trusted.
    """
    df = load_mutheta_sim(mutheta_path)
    genes = load_gene_list(gene_list_path)
    if len(genes) != len(df):
        raise ValueError(
            f"Gene list has {len(genes)} entries but sim data has {len(df)} records -- "
            "positional join is not valid."
        )
    df = df.copy()
    df["gene"] = genes
    return df


def real_post_transcriptional_regulator(
    genes: list[str] | None = None,
    connection_matrix_path: str = "data/raw/Connection_Matrix.xlsx",
) -> dict[str, str]:
    """Maps gene -> its canonical post-transcriptional regulator label, for
    genes targeted by *exactly one* of the 6 real post-transcriptional
    regulators (ambiguous genes targeted by 2+ are excluded, since there's
    no single label to test against)."""
    from data_loader import REGULATOR_INFO, load_clean

    targets, _, _ = load_clean(connection_matrix_path)
    targets = targets.reindex(REGULATOR_ORDER)

    if genes is None:
        genes = list(dict.fromkeys(load_gene_list()))  # de-duplicated, order preserved

    result = {}
    for gene in genes:
        if gene not in targets.columns:
            continue
        hits = [lab for lab in POST_TRANSCRIPTIONAL_LABELS if targets.loc[lab, gene] == 1]
        if len(hits) == 1:
            result[gene] = hits[0]
    return result


def best_theta_to_regulator_mapping(
    mutheta_path: str = DEFAULT_PATH,
    gene_list_path: str = GENE_LIST_PATH,
    connection_matrix_path: str = "data/raw/Connection_Matrix.xlsx",
) -> dict:
    """Brute-force tests every permutation of theta_0..theta_5 against the 6
    real post-transcriptional regulator labels, scoring each by per-gene
    prediction accuracy (argmax of that gene's mean Theta vector, across its
    sim rows, vs. its real single post-transcriptional regulator).

    Returns a dict with: ``n_genes_evaluated``, ``majority_baseline_acc``
    (always guessing the most common real label), ``best_perm`` (theta_i ->
    label list) and ``best_acc``, ``worst_perm``/``worst_acc``, and
    ``mean_acc_all_permutations`` (should be ~1/6 = chance).

    See module docstring: the empirical result is that ``best_acc`` is
    *below* ``majority_baseline_acc`` -- Theta's argmax does not identify a
    gene's real post-transcriptional regulator better than guessing lhp-1
    for everyone.
    """
    sim = load_mutheta_sim_with_genes(mutheta_path, gene_list_path)
    theta_cols = [f"theta_{i}" for i in range(N_THETA)]
    mean_theta_by_gene = sim.groupby("gene")[theta_cols].mean()

    real_regulator = real_post_transcriptional_regulator(
        genes=list(mean_theta_by_gene.index), connection_matrix_path=connection_matrix_path
    )
    eval_genes = list(real_regulator.keys())
    y_true = [real_regulator[g] for g in eval_genes]
    argmax_idx = mean_theta_by_gene.loc[eval_genes, theta_cols].values.argmax(axis=1)

    majority_label = pd.Series(y_true).value_counts().idxmax()
    majority_baseline_acc = sum(y == majority_label for y in y_true) / len(y_true)

    scored = []
    for perm in itertools.permutations(POST_TRANSCRIPTIONAL_LABELS):
        pred = [perm[i] for i in argmax_idx]
        acc = sum(p == t for p, t in zip(pred, y_true)) / len(y_true)
        scored.append((acc, list(perm)))
    scored.sort(key=lambda x: x[0], reverse=True)

    return {
        "n_genes_evaluated": len(eval_genes),
        "majority_baseline_acc": majority_baseline_acc,
        "majority_label": majority_label,
        "best_perm": scored[0][1],
        "best_acc": scored[0][0],
        "worst_perm": scored[-1][1],
        "worst_acc": scored[-1][0],
        "mean_acc_all_permutations": sum(s[0] for s in scored) / len(scored),
    }


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
    df = load_mutheta_sim_with_genes()
    print(f"Loaded {len(df)} records, {df['gene'].nunique()} distinct genes")
    print(df.head())
    print()
    print("Dominant theta column distribution:")
    print(df["dominant_theta"].value_counts(normalize=True).round(3))
    print()
    print("Comparison to real post-transcriptional regulator target shares:")
    print(compare_theta_dominance_to_real_network().round(3))
    print()
    print("Per-gene validation against real post-transcriptional regulator identity:")
    result = best_theta_to_regulator_mapping()
    print(f"  genes evaluated: {result['n_genes_evaluated']}")
    print(f"  majority-class baseline (always '{result['majority_label']}'): {result['majority_baseline_acc']:.3f}")
    print(f"  best of 720 permutations: {result['best_acc']:.3f}  {result['best_perm']}")
    print(f"  mean over all permutations (chance): {result['mean_acc_all_permutations']:.3f}")
    if result["best_acc"] < result["majority_baseline_acc"]:
        print("  -> Theta's argmax does NOT beat the majority baseline. Hypothesis refuted.")
