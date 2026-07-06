"""Load and clean the VTENS-fitted regulator-target connection matrix.

Source: data/raw/Connection_Matrix.xlsx, a (11 regulators x 3402 target
columns) binary matrix exported from the ensemble fit described in
Al-Omari et al. 2022 (IEEE Access 10:32510-32524).

Raw-file quirks handled here (discovered by direct inspection, not assumed
from the paper):

1. Row label whitespace: "PostReg7 NCU03363 " has a trailing space.
2. Column label whitespace: one target column is "NCU03363 " (trailing
   space), a near-duplicate of clean column "NCU03363".
3. True duplicate target columns: when openpyxl/pandas reads a sheet with
   repeated column headers, pandas mangles the later ones with ".1", ".2",
   etc. suffixes. Un-mangling (stripping a trailing ".<int>" when the base
   name is also present among the columns) recovers 3030 distinct base
   names out of 3402 raw columns -- i.e. 372 raw columns are duplicates of
   an earlier column, spread across 362 distinct genes (most duplicated
   twice; a few, including the three core-clock genes below, up to 5x).
4. Non-NCU "target" columns: three core clock genes appear with positional
   suffixes -- "wc-2_595_629", "frq_2481_2515",
   "wc-1(825-37223790-4395)_2861_2895" -- each repeated (pandas-mangled)
   5 times. These are the core clock mechanism (wc-1/wc-2/frq) taken as
   given in the paper's model, not downstream ccg targets, so they are
   excluded from the target-gene universe and returned separately.

Cleaning rule applied:
- Strip whitespace from row and column labels.
- Un-mangle pandas' duplicate-column suffixes back to base names.
- Split columns into: NCU-pattern target genes vs. core-clock-gene columns.
- Collapse duplicate NCU target columns via logical OR across each
  regulator's row (a regulator is considered to target a gene if ANY of
  that gene's duplicate columns show a 1).
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import pandas as pd

NCU_PATTERN = re.compile(r"^NCU\d{5}$")
_SUFFIX_PATTERN = re.compile(r"^(.*)\.(\d+)$")

REGULATOR_INFO = {
    "Wcc": {"gene": "WCC", "ncu": None, "class": "master"},
    "Reg1 NCU07392": {"gene": "pro-1", "ncu": "NCU07392", "class": "transcriptional"},
    "Reg2 NCU01640": {"gene": "rpn-4", "ncu": "NCU01640", "class": "transcriptional"},
    "Reg3 NCU06108": {"gene": "NCU06108", "ncu": "NCU06108", "class": "transcriptional"},
    "Reg4 NCU07155": {"gene": "NCU07155", "ncu": "NCU07155", "class": "transcriptional"},
    "PostReg5 NCU08295": {"gene": "lhp-1", "ncu": "NCU08295", "class": "post-transcriptional"},
    "PostReg6 NCU04504": {"gene": "rrp-3", "ncu": "NCU04504", "class": "post-transcriptional"},
    "PostReg7 NCU03363": {"gene": "NCU03363", "ncu": "NCU03363", "class": "post-transcriptional"},
    "PostReg8 NCU04799": {"gene": "pab-1", "ncu": "NCU04799", "class": "post-transcriptional"},
    "PostReg9 NCU00919": {"gene": "rok-1", "ncu": "NCU00919", "class": "post-transcriptional"},
    "PostReg10 NCU09349": {"gene": "has-1", "ncu": "NCU09349", "class": "post-transcriptional"},
}

CORE_CLOCK_GENE_PREFIXES = ("wc-1", "wc-2", "frq")


@dataclass
class CleaningReport:
    raw_shape: tuple[int, int]
    raw_row_labels: list[str]
    n_raw_columns: int
    n_unique_base_names: int
    n_duplicated_base_names: int
    n_columns_in_duplicate_groups: int
    duplicate_group_sizes: pd.Series
    n_core_clock_columns: int
    core_clock_base_names: list[str]
    n_target_genes_final: int
    raw_row_sums: pd.Series
    cleaned_row_sums: pd.Series

    def __str__(self) -> str:  # pragma: no cover - human-readable summary
        lines = [
            f"Raw matrix shape: {self.raw_shape[0]} regulators x {self.raw_shape[1]} columns",
            f"Raw columns: {self.n_raw_columns}",
            f"Unique base gene names after un-mangling pandas dedup suffixes: {self.n_unique_base_names}",
            f"Base names that were duplicated in the raw file: {self.n_duplicated_base_names}",
            f"Raw columns absorbed into duplicate groups: {self.n_columns_in_duplicate_groups}",
            f"Core clock gene columns excluded from target universe: {self.n_core_clock_columns} "
            f"({self.core_clock_base_names})",
            f"Final cleaned target gene count: {self.n_target_genes_final}",
            "",
            "Row sums, raw vs. cleaned (regulator -> target count):",
        ]
        cmp = pd.DataFrame(
            {"raw": self.raw_row_sums, "cleaned": self.cleaned_row_sums}
        )
        lines.append(cmp.to_string())
        return "\n".join(lines)


def _unmangle_base_names(columns: pd.Index) -> list[str]:
    """Recover pandas' pre-mangling column names.

    pandas appends ``.1``, ``.2``, ... to later occurrences of a repeated
    column label on read. We reverse that only when the un-suffixed name is
    itself present among the columns (to avoid stripping a genuine
    ``.1``-like suffix that isn't actually a mangling artifact).
    """
    cols = list(columns)
    col_set = set(cols)
    base_names = []
    for c in cols:
        m = _SUFFIX_PATTERN.match(c)
        if m and m.group(1) in col_set:
            base_names.append(m.group(1))
        else:
            base_names.append(c)
    return base_names


def load_raw_matrix(path: str) -> pd.DataFrame:
    """Load the sheet with row 0 as column headers and column 0 as the row index."""
    raw = pd.read_excel(path, sheet_name=0, header=0, index_col=0)
    raw.index = raw.index.astype(str)
    raw.columns = raw.columns.astype(str)
    return raw


def clean_matrix(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, CleaningReport]:
    """Clean the raw connection matrix.

    Returns
    -------
    targets : DataFrame (11 regulators x N distinct NCU target genes), binary,
        duplicate columns collapsed via logical OR.
    core_clock : DataFrame (11 regulators x 3 core clock genes), binary,
        duplicate columns collapsed via logical OR. Kept separate; not part
        of the ccg "target" universe.
    report : CleaningReport
    """
    raw = raw.copy()
    raw_row_labels = list(raw.index)
    raw.index = raw.index.str.strip()
    raw.columns = raw.columns.str.strip()

    raw_row_sums = raw.sum(axis=1)
    raw_row_sums.index = raw.index

    base_names = pd.Index(_unmangle_base_names(raw.columns))
    vc = pd.Series(base_names).value_counts()
    duplicate_group_sizes = vc[vc > 1].sort_values(ascending=False)

    # Group raw columns by recovered base name, collapse via logical OR.
    grouped = {}
    for base, col in zip(base_names, raw.columns):
        grouped.setdefault(base, []).append(col)

    collapsed = pd.DataFrame(
        {base: raw[cols].max(axis=1) for base, cols in grouped.items()},
        index=raw.index,
    )

    is_core_clock = collapsed.columns.str.startswith(CORE_CLOCK_GENE_PREFIXES)
    is_ncu = collapsed.columns.str.match(NCU_PATTERN)
    is_other = ~is_core_clock & ~is_ncu

    core_clock = collapsed.loc[:, is_core_clock]
    targets = collapsed.loc[:, is_ncu]

    if is_other.any():
        raise ValueError(
            f"Unclassified columns found (neither NCU pattern nor core clock gene): "
            f"{collapsed.columns[is_other].tolist()}"
        )

    report = CleaningReport(
        raw_shape=raw.shape,
        raw_row_labels=raw_row_labels,
        n_raw_columns=raw.shape[1],
        n_unique_base_names=len(vc),
        n_duplicated_base_names=len(duplicate_group_sizes),
        n_columns_in_duplicate_groups=int(duplicate_group_sizes.sum()),
        duplicate_group_sizes=duplicate_group_sizes,
        n_core_clock_columns=core_clock.shape[1],
        core_clock_base_names=list(core_clock.columns),
        n_target_genes_final=targets.shape[1],
        raw_row_sums=raw_row_sums,
        cleaned_row_sums=targets.sum(axis=1),
    )

    return targets, core_clock, report


def load_clean(path: str) -> tuple[pd.DataFrame, pd.DataFrame, CleaningReport]:
    raw = load_raw_matrix(path)
    return clean_matrix(raw)


if __name__ == "__main__":
    targets, core_clock, report = load_clean("data/raw/Connection_Matrix.xlsx")
    print(report)
    print()
    print("Target matrix shape:", targets.shape)
    print("Core clock gene matrix shape:", core_clock.shape)
