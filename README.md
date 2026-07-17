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

**Done:** `src/data_loader.py`, `src/network_builder.py`, `src/ffl_analysis.py`,
`src/null_models.py`, `src/enrichment.py`, `src/network_stats.py`, `src/viz.py`,
`src/sequence_families.py`, `src/link_analysis.py`, `src/gene_layout_3d.py`,
`src/binding_strength.py`, `src/mutheta_loader.py`, `app/streamlit_app.py`, plus the
full `tests/` suite (72 tests, passing).

**Not yet built:** exploration notebook. **Blocked on external data:** interpreting
`data/raw/MuThetaDataSim_m4.txt` (gene/sample identifiers pending from the data
source — see "Pending: simulated Mu/Theta dataset" below).

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
regulated by **at most 2** of the 11 regulators. The paper's underlying model
(**Variable Topology Ensemble Methods**, VTENS) does not fit each gene 11
independent per-regulator binding strengths, later tested one at a time. Under
VTENS's "supernet" formulation, every regulator is initially allowed to regulate
every candidate target, and the network *topology* and its *kinetic parameters*
are estimated **jointly** by MCMC — each accepted sample is one entire candidate
network consistent with the mRNA time-course data, not 11 independent draws per
gene. `Connection_Matrix.xlsx` is a thresholded snapshot of that joint fit (the
dominant regulator(s) per gene), not the full multi-way estimate. **The reported
binding strengths are therefore joint ensemble point estimates, not independent
per-parameter distributions — there is no marginal "does the ensemble exclude
zero" test to run on them**, and `data/raw/binding_strength_by_function.csv` /
the app's Binding Strengths tab exist specifically to show a real slice of these
continuous estimates without implying otherwise. That single structural fact
explains essentially every downstream mismatch: fewer multiply-regulated genes,
zero genes with 2+ *additional* regulators beyond WCC (required for an
"overlapping FFL"), and a lower total FFL count.

This is flagged rather than silently reconciled, per the project's guiding principle:
report discrepancies, don't force numbers to match.

## Novel extension: binding strength by function (`src/binding_strength.py`)

A second user-supplied file, `data/raw/binding_strength_by_function.csv`
(transcribed from the paper's own supplementary binding-strength table), gives a
real, continuous, function-level look at the joint VTENS/MCMC estimates
`Connection_Matrix.xlsx` was thresholded from: 21 KEGG/GO functional categories
(ribosome, proteasome, carbon metabolism, ...) x the same 11 regulators. The
Binding Strengths tab renders it as a heatmap (sequential blue for magnitude;
exact-zero cells render as blank surface, not pale blue, since there's a real gap
in this data between 0 and the smallest observed positive strength) plus the raw
table, with the module docstring and in-app copy both spelling out explicitly why
a `0.0000` cell means "not the estimated dominant regulator for this function" and
**not** "this regulator's ensemble distribution excluded zero" — that's not a
claim this data supports.

## Pending: simulated Mu/Theta dataset (`src/mutheta_loader.py`)

A third user-supplied file, `data/raw/MuThetaDataSim_m4.txt`, has been added and
parsed but is **not yet interpretable** — gene/sample identifiers for it are
expected from the data source in the future; until then it can't be joined to any
other file here by gene.

What's confirmed by direct inspection of the raw file (not assumed):
- The file uses bare `\r` (CR-only) line endings — splitting on `\n` collapses it
  to one line. Correctly split, it holds **2418 records**, each 16
  `index\tvalue` pairs (indices `0`-`15`, no header, no ID column).
- Indices 0-9 ("Mu"): 10 continuous values per record, ~`[0, 100]`.
- Indices 10-15 ("Theta"): 6 non-negative values that sum to **exactly 1.0 in
  every one of the 2418 records** — a probability/mixture-weight vector over 6
  categories.

That 6-category count matches this project's 6 post-transcriptional "RNA operon"
regulators (`PostReg5`.."PostReg10" in `REGULATOR_ORDER`), and one category (theta_0)
dominates the argmax far more than uniform (35.6% vs. 16.7%) — echoing lhp-1's
outsized real target share among those 6 regulators (53.9%). But
`compare_theta_dominance_to_real_network()` shows the resemblance stops there: the
relative ranking of the other 5 theta categories does **not** match the other 5
regulators' real target-count shares (e.g. theta_1's simulated dominance share is the
*smallest* of the six, while the real PostReg6/rrp-3 share is the *second largest*).
So this is flagged as a structural resemblance worth watching, not a confirmed
identity mapping — `load_mutheta_sim()` keeps the theta columns generically named
(`theta_0`..`theta_5`) by default, and only labels them as specific regulators if you
explicitly opt in (`label_theta_as_post_transcriptional=True`), which the docstring
is explicit is an assumption. A cleaned, parsed copy is cached at
`data/processed/mutheta_sim_m4.csv`.

## Novel extension: is the FFL count actually significant?

The paper's own standard of evidence for "71 FFLs is meaningful" was a single
qualitative comparison to a textbook number from a different organism (71 > 40 in
*E. coli*, Milo et al. 2002). It never asked: given *this* network's regulator target
counts, how many FFLs would you expect from a randomly rewired network with the same
degree sequence?

`src/null_models.py` answers this with a degree-preserving double-edge-swap null
model: 1,000 randomized networks (16,275 swaps each — a 5x-edges mixing budget),
holding every regulator's target count and every gene's regulator-count fixed (the
fixed WCC→regulator R structure is never touched), recomputing the FFL count each
time.

| Quantity | Value |
|---|---|
| Observed FFL count (this export) | 30 |
| Null distribution mean ± SD | 32.01 ± 4.94 |
| Z-score | −0.407 |
| Empirical p-value (upper-tail) | 0.690 |
| Empirical p-value (two-sided) | 0.775 |

**Finding: the observed FFL count is not distinguishable from a randomly rewired
network with the same degree sequence.** This is a genuinely different conclusion
from the paper's implied significance claim — but it's consistent with the earlier
finding that this export caps every gene at ≤2 regulators: with so little room for a
gene to receive both a WCC edge and a second regulator's edge, the FFL count is
mostly a mechanical consequence of the degree sequence itself, which is exactly what
the null model holds fixed. Whether the paper's *own* full continuous-ensemble network
(768 targets on lhp-1 alone, up to 11 regulators/gene) would still show this
degree-sequence-driven non-significance is an open question this export cannot answer.

## Novel extension: pairwise and triple-wise regulator co-targeting enrichment

`src/enrichment.py` tests all C(11,2) = 55 regulator pairs and C(11,3) = 165 regulator
triples for shared-target enrichment against an exact hypergeometric null (pairs) and
an exact convolved-hypergeometric null (triples — see module docstring for the exact
two-stage derivation), then applies Benjamini–Hochberg FDR correction across all
55 + 165 = 220 tests at α = .05 (full results: `data/processed/enrichment_results.csv`).

**63 / 220 tests are significant after correction — and every single one of them is a
depletion, not an enrichment.** E.g. NCU06108 + lhp-1 share only 17 target genes where
~141 would be expected by chance (q ≈ 1.9e-51); no pair or triple shows significantly
*more* co-targeting than chance. This is the same structural signature as the FFL
result above, seen from a different angle: because this export assigns each gene to
at most 2 regulators, any two regulators' target sets are pushed toward mutual
exclusivity purely as a side effect of that binarization, not because of biologically
meaningful antagonism. This is flagged as a property of the data export, not
interpreted as a biological claim about the real regulatory network.

## Novel extension: sequence-based paralog gene families (`src/sequence_families.py`)

A user-supplied file, `data/raw/uniprotkb_proteome_UP000001805_2026_07_02.fasta`
(UniProt's full *N. crassa* reference proteome, 10,255 protein sequences), enables a
question neither the paper nor any of the analysis above touches: **which
clock-network genes have close evolutionary relatives (paralogs), and how large are
those paralog families?** This is the kind of gene-family definition needed to apply
a duplication-aware framework (the companion Scruse/Arnold/Robinson
Partial-Duplication project) to this network — not possible without real sequence
data.

**Method, two stages** (no BLAST/DIAMOND/MMseqs/CD-HIT available in this
environment): (1) a cheap protein 5-mer Jaccard-similarity prefilter via an
inverted index, just to shrink ~10,255² possible pairs down to a few hundred
plausible candidates; (2) real local alignment confirmation on every candidate
(Biopython, BLOSUM62, affine gaps — the same scoring approach BLAST itself uses),
requiring both percent identity and alignment coverage to clear a threshold before
two genes count as the same family. A first pass using only the k-mer stage (no
alignment confirmation) was tried and rejected: at a strict threshold it mostly
recovered *duplicate UniProt records of the same gene* (e.g. `mtA-1`/`matA-1`,
`npl4`/`NCU02680` — identical proteins filed under two different gene symbols), not
real paralogs; those duplicate-record pairs are explicitly detected (≥97% identity,
≥85% coverage) and collapsed to one canonical gene *before* real family clustering,
so family sizes reflect genuine duplication, not annotation redundancy.

**Regulator ID resolution surfaced real discrepancies, investigated directly against
the raw file rather than assumed:**

| Regulator (paper name) | NCU ID | Resolution |
|---|---|---|
| WCC | — | Excluded — WCC is the WC-1/WC-2 protein complex, not a single gene product |
| pro-1 | NCU07392 | **Unresolved** — UniProt's own `GN=pro-1` (Q12641) is a well-established, distinct proline-biosynthesis enzyme, not a DNA-binding regulator; likely a gene-symbol collision, not verified as the same locus |
| NCU03363 | NCU03363 | **Unresolved** — no matching entry found by locus tag or any name variant tried; possibly renumbered in a genome reannotation since 2022 |
| pab-1 | NCU04799 | Resolved via alias — UniProt's `GN=pabp-1` (Q7S6N6) is the same protein |
| rok-1 | NCU00919 | Resolved via alias — UniProt's `GN=drh-16` (Q7SFC8, entry name ROK1_NEUCR) is the same protein |
| rpn-4, has-1, rrp-3, lhp-1, NCU06108, NCU07155 | — | Resolved directly |

8 of 11 regulators were confidently resolved to a sequence; 3 were not, and are
excluded from family analysis rather than guessed.

**Result:** none of the 8 resolved regulators have a detected paralog at all (each is
its own singleton "family of 1") — this analysis found no evidence that any of the
paper's clock regulators belongs to a small paralog family. Among the ~3,026 *target*
genes, real (biologically plausible) families did turn up, including:

| k | Families found | Example |
|---|---|---|
| 3 | 6 | `NCU08989`/`NCU07173`/`NCU08340` — all ADP-ribosylation factors; `NCU08340` is a clock-network target gene |
| 4 | 5 | `gh11-1`/`gh11-2`/`NCU05159`/`NCU09664` — xylanase/xylan-esterase family; `NCU05159` is a clock-network target gene |

Full gene→family assignments (9,747 genes after duplicate-record collapsing; family
size distribution `{1: 9538, 2: 72, 3: 6, 4: 5, 5: 1, 6: 1, 8: 2}`) are cached at
`data/processed/gene_families.csv` and browsable in the Streamlit app's **Gene
Families** tab, including every clock-network gene that has *any* detected paralog
(not just k=3/4).

**Known limitations, stated directly rather than glossed over:** the k-mer prefilter
can miss true but highly divergent (remote) paralogs entirely (no profile/HMM step);
single-linkage clustering of confirmed edges can over-merge distinct families through
a promiscuous shared domain; there's no e-value-style statistical significance behind
the identity/coverage thresholds, unlike real BLAST. Treat every family here as a
first-pass, biologically plausible hypothesis, not a publication-grade orthology call.

## Descriptive network statistics (`src/network_stats.py`)

Hub sizes, the singly-/multiply-regulated gene split, and a log-log power-law fit
over the combined degree distribution (11 regulator nodes + target-gene nodes),
compared to the paper's Section V-E:

| Quantity | Paper | This export | Match? |
|---|---|---|---|
| Largest hub | lhp-1, 768 | lhp-1, 739 | Close (same regulator) |
| Smallest hub | rok-1, 54 | rok-1, 54 | Exact |
| Avg. node degree | 2.23 | 2.14 | Close |
| Power-law exponent | −1.4 | −1.16 (R²=0.85) | Close, given caveat below |
| Singly-regulated genes | 2686 | 2797 | Close, not exact |
| Multiply-regulated genes | 694 | 229 | Diverges (see above) |

Interestingly, the coarse shape of the degree distribution (average node degree,
power-law exponent) tracks the paper reasonably closely even though the
multiply-regulated gene count doesn't — but the fit here is drawn from a much
poorer sample (this export's per-gene degree only takes values `{1, 2}`, vs. the
paper's up to 11), so treat the exponent as illustrative rather than a robust
independent confirmation.

## Visualization (`src/viz.py`, `src/link_analysis.py`, `src/gene_layout_3d.py`) and the Streamlit app (`app/streamlit_app.py`)

`viz.py` builds the Plotly figures — reused identically by the CLI scripts and the
app — using a fixed, colorblind-validated palette: **`REGULATOR_COLORS`**, one
distinct hue per regulator (11 total, extending the base 8-slot categorical palette
with 3 more hues chosen for CVD separation), used consistently everywhere a
regulator is drawn; plus the blue↔red diverging pair for enriched/depleted, and one
sequential blue hue for magnitude/fit lines and the binding-strength heatmap:
- A WCC-centered hub-and-spoke diagram (paper Fig. 1A layout, 2D or 3D): gray spokes
  are the fixed WCC→regulator structure; colored chords mark pairwise enrichment
  (blue) or depletion (red) significant at an adjustable q-value threshold.
  `hub_network_figure`/`hub_network_figure_3d` take `q_threshold`,
  `show_enriched`/`show_depleted`, and `focus_label` params so the app can re-filter
  and highlight interactively without recomputing any statistics (a BH-adjusted
  q-value already *is* "significant at any alpha", so re-thresholding is just
  `q_value <= new_threshold` — no re-correction needed).
- **Regulators + all 3,026 target genes**, force-directed, in two engines: a 2D graph
  via Cytoscape.js's `fcose` layout (`link_analysis.py` builds the elements/styles,
  rendered by the `st_link_analysis` component — click a node to highlight its
  neighbors, built into the component), and a true 3D layout via networkx's
  `spring_layout(dim=3)` (`gene_layout_3d.py`, precomputed/cached since it takes ~45s
  for this graph, rendered by `viz.regulator_gene_figure_3d`). Same per-regulator
  coloring in both.
- A feedforward-loop network diagram (`ffl_network_figure`) — the same FFL count from
  the FFLs & Null Model tab, laid out as the actual WCC→B→C / WCC→C loops it's
  counting.
- The log-log degree-distribution plot with the fitted exponent vs. the paper's −1.4.
- The null-model FFL-count histogram with the observed value marked.
- A binding-strength-by-function heatmap (`binding_strength_heatmap`) — sequential
  blue for magnitude, exact zeros rendered as blank surface.

`app/streamlit_app.py` is an 8-tab dashboard — **Overview**, Matrix Explorer,
Binding Strengths, FFLs & Null Model, Enrichment Tests, Network Viz, Gene Families,
Methodology & Paper Comparison — with the tabs immediately below the title (no long
blurb blocking them) and every tab written for someone who hasn't read the paper:
- **Overview** carries the explanatory blurb (what the clock network is, what the
  paper found, what this project adds, the binarized-export caveat, and the joint
  VTENS/MCMC framing behind it) plus a "how to use this app" walkthrough and a
  glossary (regulator, RNA operon, FFL, null model, Z-score, hypergeometric test,
  q-value, enriched/depleted) covering every term used elsewhere in the app.
- **Matrix Explorer**'s gene lookup is a searchable dropdown (type to filter by NCU
  ID) rather than free text, so an unfamiliar user can't mistype a gene ID into a dead
  end.
- **Binding Strengths** opens with an explicit "what this is *not*" callout: these are
  joint ensemble point estimates, not independent per-parameter distributions, and a
  zero cell is not the result of an ensemble-excludes-zero test.
- **Network Viz** has three subtabs — **Hub structure** (the 2D/3D toggle above, with
  a "Focus on a regulator" dropdown showing a paper-sourced description
  (`REGULATOR_BLURBS`, drawn from Sections V-B/V-C/V-F), hub size, FFL contribution,
  and significant relationships), **Regulators + Genes** (the force-directed 2D/3D
  graph above), and **FFL Network** (the loop diagram above).
- **Gene Families** shows the regulator ID-resolution table (resolved/aliased/
  unresolved, with reasons), the overall family-size distribution, every k=3/k=4
  family found with a "touches the clock network" marker, a browsable table of every
  clock-network gene with any detected paralog, and a searchable gene→family lookup.
- Every tab opens with an `st.expander("ℹ️ ...")` explaining what the controls do and
  what to look for, so no tab assumes prior context from another one.

It reads the cached results (`data/processed/null_model_results.npz`,
`data/processed/enrichment_results.csv`, `data/processed/gene_families.csv`,
`data/processed/gene_layout_3d.npz`) rather than recomputing them live, since the
null model alone takes ~8-9 minutes, the gene-family pipeline ~3 minutes, and the 3D
layout ~45s. Run it with:

```
streamlit run app/streamlit_app.py
```

A `.streamlit/config.toml` sets a light theme matching the same palette used in the
figures (background `#fcfcfb`, primary accent `#2a78d6`).

## Methodology: reproduced vs. novel

- **Reproduced from the paper**: R/T matrix construction (Section III, "the regulatory
  CCG proteins do not regulate each other"), the FFL-counting formula via the Hadamard
  product `(R @ T) * T` (Section V-F), the per-regulator "B-gene" breakdown via
  inner products of the WCC row of T against each other row (explicitly stated in
  Section V-F), and the function-level binding-strength table (the paper's own joint
  VTENS/MCMC point estimates, transcribed as supplied).
- **Novel extension**: a degree-preserving double-edge-swap null model giving an
  empirical Z-score/p-value for the FFL count (`null_models.py`); exact
  hypergeometric pairwise and convolved-hypergeometric triple-wise regulator
  co-targeting enrichment tests across all 55 + 165 = 220 combinations, BH-FDR
  corrected at α = .05 (`enrichment.py`). The published paper never computed a
  null-model p-value for any of these quantities — both extensions surface a
  structural property of this particular export (near-mutual-exclusive regulator
  assignment) that a purely qualitative comparison would have missed.
- **Novel extension, different data source entirely**: sequence-based paralog gene
  families from the UniProt reference proteome FASTA, using a two-stage k-mer
  prefilter + real BLOSUM62 alignment confirmation, to identify clock-network genes
  belonging to small (k=3, k=4) paralog families (`sequence_families.py`) — something
  neither the paper nor any of the network-topology analysis above could answer, since
  none of it touches actual gene sequences.

## Repo layout

```
neurospora_clock_model/
├── .streamlit/config.toml    # theme, for local `streamlit run`
├── data/
│   ├── raw/
│   │   ├── Connection_Matrix.xlsx
│   │   ├── binding_strength_by_function.csv
│   │   ├── MuThetaDataSim_m4.txt      # pending gene IDs from data source, see README
│   │   └── uniprotkb_proteome_UP000001805_2026_07_02.fasta
│   └── processed/             # cached null-model, enrichment, gene-family, 3D-layout, mutheta-sim results (app reads these)
├── src/
│   ├── data_loader.py         # done
│   ├── network_builder.py     # done
│   ├── ffl_analysis.py        # done
│   ├── null_models.py         # done
│   ├── enrichment.py          # done
│   ├── network_stats.py       # done
│   ├── viz.py                 # done
│   ├── link_analysis.py       # done
│   ├── gene_layout_3d.py      # done
│   ├── binding_strength.py    # done
│   ├── mutheta_loader.py      # done, not yet wired into app -- see README
│   └── sequence_families.py   # done
├── app/streamlit_app.py       # done
├── tests/                     # done (72 tests)
│   ├── test_ffl_count.py
│   ├── test_viz.py
│   ├── test_link_analysis.py
│   ├── test_sequence_families.py
│   └── test_mutheta_loader.py
└── notebooks/
```

## Running what exists

```
pip install -r requirements.txt
python src/data_loader.py       # cleaning report
python src/network_builder.py   # R/T construction
python src/ffl_analysis.py      # FFL count + validation vs. paper
python src/network_stats.py     # degree distribution, hub sizes, power-law fit
python src/null_models.py       # degree-preserving null model (~8-9 min, 1000 iterations); caches to data/processed/
python src/enrichment.py        # pairwise/triple enrichment + BH-FDR (~45s); caches to data/processed/
python src/sequence_families.py # sequence-based paralog families (~3 min); caches to data/processed/
python src/gene_layout_3d.py    # 3D force-directed regulator+gene layout (~45s); caches to data/processed/
python src/binding_strength.py  # binding-strength-by-function table, printed to stdout (no cache needed)
python -m pytest tests/ -v      # test suite (~50s)
streamlit run app/streamlit_app.py   # interactive dashboard (reads the cached results above)
```
