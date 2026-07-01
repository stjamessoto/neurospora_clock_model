# Neurospora Clock Network: Subnetwork Enrichment Model

Statistical extension of the fitted clock regulatory network from:

> Al-Omari, A. M., Griffith, J., Scruse, A., Robinson, R. W., Schüttler, H.-B., &
> Arnold, J. (2022). Ensemble Methods for Identifying RNA Operons and Regulons in the
> Clock Network of *Neurospora Crassa*. *IEEE Access*, 10, 32510-32524.
> https://doi.org/10.1109/ACCESS.2022.3160481

That paper used Variable Topology Ensemble Methods (VTENS) — MCMC fit against 60,759
mRNA time points — to reconstruct a genome-scale clock network of 3,380 genes
regulated by 12 hubs: the master regulator WCC, 5 transcriptional regulators, and 6
post-transcriptional "RNA operon" regulators (`data/raw/Connection_Matrix.xlsx` is the
exported target-assignment matrix from that fit).

This project treats that matrix as real input data and adds formal statistical testing
(null models, Z-scores, Benjamini-Hochberg FDR correction) that the original paper did
not perform — it only did qualitative size comparisons (e.g. "71 FFLs > 40 in E.
coli").

## Status

**Done (data audit through FFL validation):** `src/data_loader.py`,
`src/network_builder.py`, `src/ffl_analysis.py`.

**Not yet built:** degree-preserving null model (`null_models.py`), pairwise/triple
regulator enrichment testing with BH-FDR (`enrichment.py`), full descriptive network
statistics (`network_stats.py`), visualization (`viz.py`), Streamlit app, test suite.

## Data cleaning findings (`src/data_loader.py`)

The raw sheet is (11 regulators x 3402 target columns). Cleaning steps and what they
found, run directly against the file rather than assumed:

- **Row label whitespace**: `"PostReg7 NCU03363 "` has a trailing space. Stripped.
- **Column label whitespace**: one target column is `"NCU03363 "` (trailing space), a
  near-duplicate of the clean column `"NCU03363"`.
- **True duplicate target columns**: pandas mangles repeated column headers on read by
  appending `.1`, `.2`, ... Un-mangling (stripping a trailing `.<int>` when the
  un-suffixed name is also present) recovers **3029 distinct base names out of 3402
  raw columns** — 363 base names are duplicated, absorbing 736 raw columns. Most
  duplicated genes appear twice; the three core clock genes appear 5x each.
- **Non-NCU "target" columns**: three core clock genes appear with positional
  suffixes — `wc-2_595_629`, `frq_2481_2515`, `wc-1(825-37223790-4395)_2861_2895` —
  each duplicated 5x by the mangling above. These are the core clock mechanism
  (wc-1/wc-2/frq) taken as given in the paper's model, **not** downstream
  clock-controlled genes, so they're excluded from the target-gene universe and kept
  separately (`core_clock` DataFrame).
- **Cleaning rule applied**: strip whitespace, un-mangle duplicate columns, collapse
  duplicates via logical OR per regulator, split out core-clock-gene columns.
- **Result**: final cleaned target matrix is (11, 3026).

## Validation against the paper (reproduced vs. diverges)

| Quantity | Paper (Fig. 9 / Sec. V-E) | This export, cleaned | Match? |
|---|---|---|---|
| rok-1 (NCU00919) target count | 54 | 54 | Exact |
| lhp-1 (NCU08295) target count | 768 | 739 | Close, not exact |
| Singly-regulated genes | 2686 | 2797 | Close, not exact |
| Multiply-regulated genes | 694 | 229 | **Diverges** |
| Total clock genes | 3380 | 3026 | Close, not exact |
| Max regulators on one gene | up to 11 (implied) | **2** | **Diverges** |
| Total FFLs (`RT ∘ T` Hadamard product) | 71 | 30 | Diverges |
| Genes with overlapping FFLs (>1 B-regulator) | 5 | 0 | Diverges |
| Top FFL-contributing regulator | lhp-1 (39 combined w/ NCU06108) | lhp-1 (9, single largest) | Same qualitative pattern, different magnitude |

**Root cause of the divergence (not a bug):** every gene in this exported matrix is
regulated by **at most 2** of the 11 regulators. The paper's model assigns each target
gene a continuous *binding strength* per regulator (Table 1; strengths across
regulators sum to 1 per gene, with Monte Carlo sampling over the full ensemble). This
`Connection_Matrix.xlsx` export looks like a thresholded/binarized snapshot of that
ensemble (e.g. top-1 or top-2 dominant regulator per gene) rather than the full
multi-way soft assignment. That single structural fact explains essentially every
downstream mismatch: fewer multiply-regulated genes, zero genes with 2+ *additional*
regulators beyond WCC (required for an "overlapping FFL"), and a lower total FFL count.

This is flagged rather than silently reconciled, per the project's guiding principle:
report discrepancies, don't force numbers to match.

## Methodology: reproduced vs. novel

- **Reproduced from the paper**: R/T matrix construction (Section III, "the regulatory
  CCG proteins do not regulate each other"), the FFL-counting formula via the Hadamard
  product `(R @ T) * T` (Section V-F), and the per-regulator "B-gene" breakdown via
  inner products of the WCC row of T against each other row (explicitly stated in
  Section V-F).
- **Novel extension (planned)**: a degree-preserving double-edge-swap null model to
  get an empirical p-value for the FFL count; hypergeometric pairwise/triple-wise
  regulator co-targeting enrichment tests across all 55 + 165 = 220 combinations,
  BH-FDR corrected at α = .05. The published paper never computed a null-model
  p-value for any of these quantities.

## Repo layout

```
neurospora_clock_model/
├── data/raw/Connection_Matrix.xlsx
├── src/
│   ├── data_loader.py       # done
│   ├── network_builder.py   # done
│   ├── ffl_analysis.py      # done
│   ├── null_models.py       # not yet built
│   ├── enrichment.py        # not yet built
│   ├── network_stats.py     # not yet built
│   └── viz.py                # not yet built
├── app/streamlit_app.py      # not yet built
├── tests/                     # not yet built
└── notebooks/
```

## Running what exists

```
pip install -r requirements.txt
python src/data_loader.py       # cleaning report
python src/network_builder.py   # R/T construction
python src/ffl_analysis.py      # FFL count + validation vs. paper
```
