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
full `tests/` suite (80 tests, passing).

**Not yet built:** exploration notebook. **Open question:** which of this project's
10 non-WCC regulators each of `MuThetaDataSim_m4.txt`'s `mu_0`..`mu_9` columns
positionally corresponds to, and what the Theta block represents — gene identity is
confirmed (`ALLNCU.txt`, one-to-one, same record count as the sim file) and Mu's
general nature is now confirmed too (conditional ensemble averages, not simple
means — see below), but the leading hypothesis for Theta (per-gene
post-transcriptional-regulator assignment) has been tested and refuted — see
"Simulated Mu/Theta dataset" below.

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
dominant regulator(s) per gene), not the full multi-way estimate. Since every gene
tops out at 2 regulators, no gene can ever be a shared target of 3 or more
regulators at once — so every 3-way and 4-way overlap in this export is
*structurally* zero, not just statistically depleted (see the enrichment section
below). **The reported
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

## Simulated Mu/Theta dataset: gene identity confirmed, regulator hypothesis refuted (`src/mutheta_loader.py`)

Two more user-supplied files: `data/raw/MuThetaDataSim_m4.txt` (a simulated dataset)
and `data/raw/ALLNCU.txt` (its gene-identifier list, supplied afterward, same row
order). The two files line up record-for-record: `MuThetaDataSim_m4.txt` parses to
exactly 2418 records and `ALLNCU.txt` has exactly 2418 lines — checked directly
(`len(_read_records(...)) == len(load_gene_list(...))`), not assumed from the
file-format writeup alone. Row order is a confirmed one-to-one correspondence (per
the data source, not just inferred from the matching counts): this is the indexing
used throughout the analysis to associate each record's Mu/Theta values with its
gene.

**Format, confirmed by direct inspection of the raw file:**
- `MuThetaDataSim_m4.txt` uses bare `\r` (CR-only) line endings — splitting on `\n`
  collapses it to one line. Correctly split, it holds **2418 records**, each 16
  `index\tvalue` pairs (indices `0`-`15`, no header, no ID column of its own).
- Indices 0-9 ("Mu"): 10 continuous values per record, ~`[0, 100]`. **Confirmed:**
  these are *conditional* ensemble averages, not simple means, over an ensemble of
  2,000 MCMC runs — for each ensemble member, the regulator with the largest
  inferred binding strength for that gene is identified as dominant, and a Mu
  column's value is the average binding strength of its regulator *only* over the
  ensemble members where that regulator was dominant for the gene (not an average
  over all 2,000 members). Still open: which of this project's 10 non-WCC
  regulators each of `mu_0`..`mu_9` positionally corresponds to.
- Indices 10-15 ("Theta"): 6 non-negative values that sum to **exactly 1.0 in
  every one of the 2418 records** — a probability/mixture-weight vector over 6
  categories. The conditional-averaging confirmation above is specific to Mu;
  Theta's identity is separately tested (and refuted, for the one hypothesis
  tried) below.
- `ALLNCU.txt` is 2418 lines, 2216 distinct NCU gene IDs (202 genes repeat, each
  repeat getting its own sim row/Mu/Theta values).

**Confirmed: this is this project's gene universe.** Every one of the 2216 distinct
genes in `ALLNCU.txt` is a real target-gene column in the cleaned
`Connection_Matrix.xlsx` matrix — 100% overlap. Not a coincidental resemblance; this
dataset shares its gene space with the rest of the project.

**Tested and refuted: Theta is not a per-gene post-transcriptional-regulator
posterior.** The initial hypothesis (motivated by the matching count of 6 categories,
and theta_0 dominating the argmax similarly to how lhp-1 dominates real target
counts) was that the 6 Theta values encode each gene's soft assignment across the 6
real post-transcriptional "RNA operon" regulators (`PostReg5`.."PostReg10"). With
gene identity now available, `best_theta_to_regulator_mapping()` tests this directly:
for the 887 sim genes targeted by exactly one real post-transcriptional regulator, it
brute-forces all 6! = 720 permutations of "theta column → regulator label" and scores
per-gene prediction accuracy (argmax of that gene's mean Theta vs. its real
regulator).

| Quantity | Value |
|---|---|
| Genes evaluated | 887 |
| Majority-class baseline (always guess lhp-1) | 0.564 |
| **Best of 720 permutations** | **0.287** |
| Mean over all permutations (chance) | 0.167 |

**The best achievable accuracy (0.287) is below the trivial majority-guess baseline
(0.564)** — Theta's argmax is a *worse* predictor of a gene's real post-transcriptional
regulator than simply always guessing lhp-1. The category-count match and the
single-dominant-category shape that motivated the hypothesis were apparently
coincidental. `load_mutheta_sim()` therefore keeps the theta columns generically
named (`theta_0`..`theta_5`); the opt-in `label_theta_as_post_transcriptional=True`
labeling is kept available but its docstring is explicit that it's an unconfirmed,
tested-and-rejected assumption, not a default to build further analysis on.

What Mu and Theta actually represent remains an open question — flagged rather than
forced, per this project's convention. A cleaned, gene-joined copy is cached at
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

## Novel extension: pairwise, triple-wise, and quartet-wise regulator co-targeting enrichment

`src/enrichment.py` tests all C(11,2) = 55 regulator pairs, C(11,3) = 165 regulator
triples, and C(11,4) = 330 regulator quartets for shared-target enrichment against an
exact hypergeometric null (pairs) and an exact convolved-hypergeometric null (triples
and quartets — see module docstring for the sequential-convolution derivation, which
generalizes to any group size), then applies Benjamini–Hochberg FDR correction across
all 55 + 165 + 330 = 550 tests at α = .05 (full results:
`data/processed/enrichment_results.csv`).

**60 / 550 tests are significant after correction (41 pairs, 19 triples, 0
quartets) — and every single one of them is a depletion, not an enrichment.** E.g.
NCU06108 + lhp-1 share only 17 target genes where ~141 would be expected by chance
(q ≈ 4.6e-51); no pair, triple, or quartet shows significantly *more* co-targeting
than chance.

**Every triple and every quartet — not just the significant ones — has an observed
overlap of exactly zero, with no exceptions across all 495 of them.** This isn't a
statistical pattern that happens to hold; it's a direct, mechanical consequence of
the "at most 2 regulators per gene" finding above — a gene can only ever belong to
the intersection of 3+ regulators' target sets if it's targeted by all of them
simultaneously, which this export's binarization makes impossible by construction.
The quartet tests exist mainly to make this structural ceiling explicit and
quantify how far short of chance-level co-targeting it falls (e.g. the closest
quartet to significance, pro-1 + NCU06108 + NCU07155 + lhp-1, still only reaches
p ≈ 0.12, q ≈ 0.71 — not close, because ~2.8 shared targets would already have been
expected by chance for that foursome). At the pairwise level, some pairs *do* share
targets (45 of 55 pairs have nonzero overlap) — it's specifically going from 2 to 3+
simultaneous regulators that the export's structure forbids outright. This is
flagged as a property of the data export, not interpreted as a biological claim
about the real regulatory network.

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

`app/streamlit_app.py` is a 10-tab dashboard — **Overview**, Matrix Explorer,
Binding Strengths, **Mu/Theta Sim**, FFLs & Null Model, Enrichment Tests, Network
Viz, Gene Families, Methodology & Paper Comparison, **Full Paper** — with the tabs
immediately below the title and every tab written for someone who hasn't read the
paper. Every
explanatory block throughout the app lives inside a collapsed-by-default
`st.expander("ℹ️ ...")` rather than static text, so a tab's default view is metrics,
charts, and tables — the "why"/"how to read this" prose is a click away, not
scrolled past:
- **Overview** is now a short intro + 3 metrics up front, with "what this project
  adds," "one important caveat" (the binarized-export/joint-VTENS-MCMC framing), and
  the full "how to use this app" walkthrough behind expanders — plus a glossary
  (regulator, RNA operon, FFL, null model, Z-score, hypergeometric test, q-value,
  enriched/depleted) covering every term used elsewhere in the app.
- **Matrix Explorer**'s gene lookup is a searchable dropdown (type to filter by NCU
  ID) rather than free text, so an unfamiliar user can't mistype a gene ID into a dead
  end.
- **Binding Strengths** opens with an explicit "what this is *not*" callout: these are
  joint ensemble point estimates, not independent per-parameter distributions, and a
  zero cell is not the result of an ensemble-excludes-zero test.
- **Mu/Theta Sim** surfaces `src/mutheta_loader.py`'s findings directly: record/gene
  counts, the best-vs-majority-baseline accuracy metrics that refute the
  regulator-assignment hypothesis, a simulated-vs-real dominant-share bar chart, a
  per-gene Mu/Theta lookup dropdown, and the full joined dataset.
- **Network Viz** has three subtabs — **Hub structure** (the 2D/3D toggle above, with
  a "Focus on a regulator" dropdown showing a paper-sourced description
  (`REGULATOR_BLURBS`, drawn from Sections V-B/V-C/V-F), hub size, FFL contribution,
  and significant relationships), **Regulators + Genes** (the force-directed 2D/3D
  graph above), and **FFL Network** (the loop diagram above).
- **Gene Families** shows the regulator ID-resolution table (resolved/aliased/
  unresolved, with reasons), the overall family-size distribution, every k=3/k=4
  family found with a "touches the clock network" marker, a browsable table of every
  clock-network gene with any detected paralog, and a searchable gene→family lookup.
- **Methodology & Paper Comparison** keeps the validation table and citation visible
  (they're compact and load-bearing) but tucks the file-by-file reproduced-vs-novel
  breakdown and the "why the numbers don't match" writeup behind expanders; a button
  at the bottom jumps to the **Full Paper** tab.
- **Full Paper** is `docs/full_methodology_paper.md` (~7,000 words): a complete,
  plain-language explainer of the biology, every method, and every result in this
  project — written so someone with no statistics or molecular-biology background can
  read it start to finish, with every technical term defined before it's used. The tab
  itself shows just the document's "Executive Summary" section (extracted from the
  markdown source, not duplicated by hand) plus three buttons: download as PDF,
  download as Word (`.docx`), and **open the PDF directly in a new browser tab** —
  implemented client-side (no server route or mirrored file needed): the PDF's bytes
  are base64-embedded in the page, and a button click decodes them into a `Blob` and
  opens it via `window.open(window.URL.createObjectURL(blob))`. Two non-obvious
  browser quirks had to be worked around to get there, both confirmed by direct
  testing rather than assumed: (1) a plain `data:` URI on an `<a target="_blank">`
  is silently blocked by Chromium's top-level-navigation policy, so the click has to
  build a `blob:` URL instead, which isn't subject to that block; (2) the button has
  to live inside `st.components.v1.html` (a real iframe), not
  `st.markdown(unsafe_allow_html=True)`, because Streamlit's markdown path strips
  inline `onclick` attributes even with that flag set — and once inside an inline
  `onclick=""` attribute specifically, the browser's legacy `with(document){
  with(form){ with(element){...} } }` wrapping around event-handler attributes means
  a bare `URL` reference resolves to `document.URL` (a string) instead of the global
  `URL` constructor, so the handler has to say `window.URL.createObjectURL(...)`
  explicitly. Both PDF/Word files are generated via `src/generate_paper.py` (uses
  `pandoc` + a LaTeX engine, both external system dependencies — see that script's
  docstring).

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
  hypergeometric pairwise and convolved-hypergeometric triple- and quartet-wise
  regulator co-targeting enrichment tests across all 55 + 165 + 330 = 550
  combinations, BH-FDR corrected at α = .05 (`enrichment.py`). The published paper
  never computed a null-model p-value for any of these quantities — both
  extensions surface a structural property of this particular export
  (near-mutual-exclusive regulator assignment, and a hard ceiling of zero observed
  overlap for any 3+ regulators at once) that a purely qualitative comparison would
  have missed.
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
│   │   ├── MuThetaDataSim_m4.txt      # simulated Mu/Theta dataset -- see README
│   │   ├── ALLNCU.txt                 # gene IDs for MuThetaDataSim_m4.txt, same row order
│   │   └── uniprotkb_proteome_UP000001805_2026_07_02.fasta
│   └── processed/             # cached null-model, enrichment, gene-family, 3D-layout, mutheta-sim results (app reads these)
├── docs/
│   ├── full_methodology_paper.md    # source of the "Full Paper" tab -- plain-language write-up
│   ├── full_methodology_paper.pdf   # generated via src/generate_paper.py, downloadable in-app
│   └── full_methodology_paper.docx  # generated via src/generate_paper.py, downloadable in-app
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
│   ├── mutheta_loader.py      # done, wired into the app's Mu/Theta Sim tab
│   ├── sequence_families.py   # done
│   └── generate_paper.py      # done, regenerates docs/full_methodology_paper.{pdf,docx} via pandoc
├── app/streamlit_app.py       # done
├── tests/                     # done (80 tests)
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
python src/enrichment.py        # pairwise/triple/quartet enrichment + BH-FDR (~1 min); caches to data/processed/
python src/sequence_families.py # sequence-based paralog families (~3 min); caches to data/processed/
python src/gene_layout_3d.py    # 3D force-directed regulator+gene layout (~45s); caches to data/processed/
python src/binding_strength.py  # binding-strength-by-function table, printed to stdout (no cache needed)
python src/generate_paper.py    # regenerates docs/full_methodology_paper.{pdf,docx} via pandoc (system dependency)
python -m pytest tests/ -v      # test suite (~50s)
streamlit run app/streamlit_app.py   # interactive dashboard (reads the cached results above)
```
