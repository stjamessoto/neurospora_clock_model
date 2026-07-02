"""Interactive dashboard for the Neurospora crassa clock network project.

Reads the fitted regulator -> target matrix from Al-Omari, Griffith, Scruse,
Robinson, Schuttler & Arnold (2022), IEEE Access 10:32510-32524, and adds
statistical tests (a degree-preserving null model for the FFL count;
BH-FDR-corrected pairwise/triple regulator co-targeting enrichment) that the
paper itself did not perform.

Run with: streamlit run app/streamlit_app.py
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import streamlit as st

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from data_loader import REGULATOR_INFO, load_clean  # noqa: E402
from network_builder import build_network  # noqa: E402
from ffl_analysis import count_ffls, PAPER_REPORTED_FFL_COUNT  # noqa: E402
from network_stats import (  # noqa: E402
    compute_network_stats,
    PAPER_AVG_NODE_DEGREE,
    PAPER_MULTIPLY_REGULATED,
    PAPER_POWER_LAW_EXPONENT,
    PAPER_SINGLY_REGULATED,
    PAPER_TOTAL_GENES,
)
from null_models import DEFAULT_CACHE_PATH as NULL_CACHE_PATH, load_result as load_null_result  # noqa: E402
from enrichment import DEFAULT_OUTPUT_PATH as ENRICHMENT_CACHE_PATH  # noqa: E402
from viz import degree_distribution_figure, hub_network_figure, null_distribution_figure  # noqa: E402
from sequence_families import (  # noqa: E402
    DEFAULT_OUTPUT_PATH as FAMILIES_CACHE_PATH,
    EXCLUDED_REGULATORS,
    REGULATOR_ALIASES,
    UNRESOLVED_REGULATORS,
    resolve_regulator_gene_names,
)

DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "Connection_Matrix.xlsx")

PAPER_CITATION = (
    "Al-Omari, A. M., Griffith, J., Scruse, A., Robinson, R. W., Schuttler, H.-B., "
    "& Arnold, J. (2022). Ensemble Methods for Identifying RNA Operons and Regulons "
    "in the Clock Network of *Neurospora Crassa*. *IEEE Access*, 10, 32510-32524. "
    "https://doi.org/10.1109/ACCESS.2022.3160481"
)

# Short, paper-grounded blurbs for each regulator, shown when a user focuses on
# it in the Network Viz tab. Facts drawn from the paper's Sections V-B/V-C/V-F.
REGULATOR_BLURBS = {
    "Wcc": (
        "The White Collar Complex (WC-1/WC-2) — the master regulator at the top of "
        "the hierarchy. It directly senses light and drives the core circadian "
        "oscillator (*frq*/*wc-1*/*wc-2*); in this network it activates all 11 "
        "downstream regulators with a fixed binding constant (μ₀ = 1)."
    ),
    "Reg1 NCU07392": (
        "**pro-1** — a transcriptional regulator with a target profile very similar "
        "to rpn-4's, but additionally enriched for ribosome, ribosome biogenesis, "
        "meiosis, and recombination genes."
    ),
    "Reg2 NCU01640": (
        "**rpn-4** — a transcriptional regulator enriched for carbon metabolism, "
        "amino acid metabolism, glycerophospholipid metabolism, and oxidative "
        "phosphorylation genes. Its RNA profile shows a strong positive response "
        "to light."
    ),
    "Reg3 NCU06108": (
        "A transcriptional regulator with a strong *negative* response to light. "
        "Along with lhp-1, it shares many ribosome/ribosome-biogenesis targets — "
        "together the two account for 39 of the paper's 71 reported feedforward "
        "loops, more than any other regulator pair."
    ),
    "Reg4 NCU07155": (
        "A transcriptional regulator the paper associates with regulation of "
        "nitrogen and sulfur metabolism."
    ),
    "PostReg5 NCU08295": (
        "**lhp-1** — the paper's single largest regulatory hub (768 targets), "
        "larger than the master regulator WCC itself. A post-transcriptional "
        "RNA-binding protein (first characterized as an RNA-operon regulator in "
        "*S. cerevisiae*) controlling mRNA *stability* for genes in the ribosome, "
        "ribosome biogenesis, proteasome, protein export, and amino acid "
        "biosynthesis; also implicated in tRNA processing."
    ),
    "PostReg6 NCU04504": (
        "**rrp-3** — an RNA-operon (post-transcriptional) regulator with a target "
        "profile similar to has-1's: ribosome, ribosome biogenesis, RNA transport, "
        "carbon metabolism, and amino acid biosynthesis."
    ),
    "PostReg7 NCU03363": (
        "A post-transcriptional RNA-operon regulator — one of the six RNA operons "
        "hypothesized in the paper."
    ),
    "PostReg8 NCU04799": (
        "**pab-1** — a poly-A binding protein and RNA-operon regulator. The only "
        "regulator in the network uniquely targeting genes homologous to those "
        "involved in synthesis and degradation of ketone bodies."
    ),
    "PostReg9 NCU00919": (
        "**rok-1** — the paper's *smallest* regulatory hub (54 targets, exactly "
        "matched in this export). A post-transcriptional RNA-operon regulator."
    ),
    "PostReg10 NCU09349": (
        "**has-1** — a DEAD-box RNA helicase and RNA-operon regulator. Its "
        "*S. cerevisiae* homolog HAS1 is required for both 40S and 60S ribosomal "
        "subunit biogenesis; its targets here include ribosome biogenesis, "
        "meiosis, RNA transport, and the spliceosome."
    ),
}

st.set_page_config(page_title="Neurospora Clock Network", layout="wide", page_icon="🕐")


@st.cache_resource
def _load_network():
    return build_network(DATA_PATH)


@st.cache_resource
def _load_cleaning_report():
    _, _, report = load_clean(DATA_PATH)
    return report


@st.cache_data
def _load_enrichment():
    if not os.path.exists(ENRICHMENT_CACHE_PATH):
        return None
    return pd.read_csv(ENRICHMENT_CACHE_PATH)


@st.cache_resource
def _load_null_model():
    return load_null_result(NULL_CACHE_PATH)


@st.cache_resource
def _load_ffl_result(_net):
    return count_ffls(_net)


@st.cache_resource
def _load_network_stats(_net):
    return compute_network_stats(_net)


@st.cache_data
def _load_gene_families():
    if not os.path.exists(FAMILIES_CACHE_PATH):
        return None
    return pd.read_csv(FAMILIES_CACHE_PATH)


net = _load_network()
cleaning_report = _load_cleaning_report()
enrichment_df = _load_enrichment()
null_result = _load_null_model()
ffl_result = _load_ffl_result(net)
stats_result = _load_network_stats(net)
families_df = _load_gene_families()
regulator_gene_names = resolve_regulator_gene_names(REGULATOR_INFO)

# ---------------------------------------------------------------------------
# Header (kept short so the tabs are visible without scrolling)
# ---------------------------------------------------------------------------

st.title("🕐 The Neurospora crassa Clock Network")
st.caption(
    "A statistical extension of the fitted RNA-operon / transcriptional-regulon "
    "network from Al-Omari et al. (2022), *IEEE Access*. New here? Start with the "
    "**Overview** tab."
)

tab_overview, tab_matrix, tab_ffl, tab_enrich, tab_viz, tab_families, tab_method = st.tabs(
    [
        "🏠 Overview",
        "📋 Matrix Explorer",
        "🔁 FFLs & Null Model",
        "🧪 Enrichment Tests",
        "🕸️ Network Viz",
        "🧬 Gene Families",
        "📖 Methodology & Paper Comparison",
    ]
)

# ---------------------------------------------------------------------------
# Tab: Overview
# ---------------------------------------------------------------------------

with tab_overview:
    st.markdown(
        f"""
**What this is.** *Neurospora crassa* (bread mold) is the model organism used to work
out how the circadian clock operates at the molecular level. Its clock is driven by a
core feedback loop (the *frq* / *wc-1* / *wc-2* genes), whose protein product **WCC**
(White Collar Complex) acts as the master regulator broadcasting time-of-day and
light information outward to the rest of the genome.

The paper behind this project — *{PAPER_CITATION}* — used a computationally
expensive Monte Carlo method (**Variable Topology Ensemble Methods**, fit against
60,759 mRNA measurements on GPUs over a 61-day run) to reconstruct which of
**3,380 clock-controlled genes** are regulated by which of **11 regulators**
downstream of WCC — 5 classic transcription factors (green, in the Network Viz tab)
and 6 **RNA operons**: post-transcriptional regulators (orange) that control mRNA
*stability* rather than transcription, per Jack Keene's "RNA operon" hypothesis. The
paper's central, striking finding was that the largest single regulatory hub in the
entire network wasn't a transcription factor — it was **lhp-1**, an RNA-binding
protein, with more target genes than the master regulator WCC itself.

**What this project adds.** The fitted network is treated here as real input data
(`data/raw/Connection_Matrix.xlsx`). The paper reported network-motif counts (e.g.
"71 feedforward loops, versus 40 in *E. coli*") as a **qualitative** size comparison.
This project asks the question the paper didn't: *compared to what, exactly?* It adds
a formal **degree-preserving null model** (is the feedforward-loop count more than
you'd expect from a randomly rewired network with the same regulator target counts?)
and **Benjamini-Hochberg FDR-corrected hypergeometric tests** for every pairwise and
triple-wise regulator combination (do any two or three regulators share more — or
fewer — target genes than chance?).

**One important caveat, upfront.** The exported connection matrix this app reads
turns out to cap every target gene at **at most 2** regulators, whereas the paper's
underlying model allows a gene to be governed by a continuous *mixture* of binding
strengths across up to all 11 regulators. That single structural difference — a
binarized/thresholded snapshot rather than the full ensemble — explains essentially
every place the numbers in this app diverge from the paper's published figures. Every
such divergence is reported explicitly rather than adjusted away; see the
**Methodology & Paper Comparison** tab for the full accounting.
"""
    )

    st.markdown("**How to use this app:**")
    st.markdown(
        """
1. **📋 Matrix Explorer** — browse the cleaned data itself: which genes each
   regulator targets, and a lookup for any individual gene.
2. **🔁 FFLs & Null Model** — a specific 3-gene network pattern (the feedforward
   loop) and whether it shows up more than random chance would predict.
3. **🧪 Enrichment Tests** — do pairs/triples of regulators share more (or fewer)
   target genes than chance? A filterable results table.
4. **🕸️ Network Viz** — an interactive diagram of the whole network; click through
   regulators to see their role and relationships.
5. **📖 Methodology & Paper Comparison** — the technical detail: every number here
   side-by-side with the paper's, and why they differ.
"""
    )

    with st.expander("📚 Glossary — key terms used throughout this app"):
        st.markdown(
            """
- **Regulator** — one of the 11 genes/proteins (plus WCC) that control other genes
  in this network. Split into two kinds:
  - **Transcriptional regulator** — a classic transcription factor; turns target
    genes' transcription on/off.
  - **RNA operon** (post-transcriptional regulator) — an RNA-binding protein that
    doesn't touch transcription, but instead controls how stable/long-lived a
    target gene's mRNA is, effectively co-regulating a whole "operon" of otherwise
    unlinked genes at the RNA level.
- **Target gene / clock-controlled gene (ccg)** — a gene downstream of a regulator,
  shown by the data to be circadian- and/or light-regulated.
- **Hub** — a regulator and all of its target genes, treated as a unit; hub *size*
  = number of target genes.
- **Feedforward loop (FFL)** — a 3-node pattern: WCC regulates a regulator, and
  both WCC and that regulator regulate the same target gene. Thought to introduce
  a delay or reinforcement effect in gene expression timing.
- **Degree** — the number of connections a node (regulator or gene) has. A
  regulator's degree = its hub size; a gene's degree = how many regulators target it.
- **Null model** — a randomized version of the same network (same regulator target
  counts, same gene regulator-counts) used as a "chance" baseline for comparison.
- **Z-score** — how many standard deviations the real (observed) value is from the
  null model's average. Roughly, |Z| > ~2 is notable; near 0 means "no different
  from chance."
- **Hypergeometric test** — the standard statistical test for "do these two sets
  share more/fewer members than you'd expect by chance," given their sizes and the
  total population.
- **p-value** — the probability of seeing an overlap this extreme (or more) by pure
  chance, if there's really no relationship.
- **q-value (BH-FDR corrected)** — a p-value adjusted for running many tests at
  once (220, here), so "5% of my *significant* results are false positives" stays
  true even after testing every regulator pair and triple. Always use q-value, not
  raw p-value, to decide significance when many tests are run together.
- **Enriched vs. depleted** — enriched = two regulators share *more* targets than
  chance predicts; depleted = *fewer* than chance predicts.
"""
        )

# ---------------------------------------------------------------------------
# Tab: Matrix explorer
# ---------------------------------------------------------------------------

with tab_matrix:
    st.subheader("Raw connection matrix, cleaned")
    with st.expander("ℹ️ How to read this tab", expanded=True):
        st.markdown(
            """
This tab shows the data behind everything else in the app: which of the 11
regulators targets which of the ~3,000 clock-controlled genes.

- The **raw vs. cleaned** table below shows how many targets each regulator has,
  before and after fixing whitespace/duplicate-column issues in the original file
  (see `src/data_loader.py`).
- Use the **gene lookup dropdown** to pick any target gene and see which
  regulator(s) control it.
- The **full matrix** at the bottom is the actual 11 x ~3,000 binary table (1 =
  "this regulator targets this gene") that every other tab is computed from.
"""
        )
    st.markdown(
        "The raw Excel export is (11 regulators x 3402 columns) with whitespace, "
        "pandas-mangled duplicate columns, and 3 core-clock-gene columns mixed in "
        "with the target genes. See `src/data_loader.py` for the full cleaning rule."
    )

    col1, col2, col3 = st.columns(3)
    col1.metric("Raw columns", cleaning_report.n_raw_columns)
    col2.metric("Cleaned target genes", cleaning_report.n_target_genes_final)
    col3.metric("Core clock genes excluded", cleaning_report.n_core_clock_columns)

    st.markdown("**Regulator target counts (row sums), raw vs. cleaned:**")
    cmp = pd.DataFrame({"raw": cleaning_report.raw_row_sums, "cleaned": cleaning_report.cleaned_row_sums})
    cmp.index = [f"{lab}  —  {REGULATOR_INFO.get(lab.strip(), {}).get('gene', '?')}" for lab in cmp.index]
    st.dataframe(cmp, width="stretch")

    st.divider()
    st.markdown("**Look up a target gene**")
    st.caption(
        "Pick a gene from the dropdown (type to search/filter by NCU ID) to see which "
        "regulator(s) target it."
    )
    gene_options = sorted(net.T.columns)
    gene_query = st.selectbox(
        "Target gene (NCU ID)",
        options=gene_options,
        index=None,
        placeholder="Start typing an NCU ID, e.g. NCU08120…",
        help="This list is every one of the cleaned target genes. Type to filter it.",
    )
    if gene_query:
        regs = net.T.index[net.T[gene_query].astype(bool)]
        reg_names = [f"{REGULATOR_INFO[r]['gene']} ({r})" for r in regs]
        st.success(f"{gene_query} is targeted by {len(reg_names)} regulator(s): {', '.join(reg_names)}")

    st.divider()
    st.markdown("**Full cleaned T matrix (regulator x target gene, binary)**")
    st.caption("Rows = regulators, columns = target genes, 1 = regulated. Scroll to explore.")
    st.dataframe(net.T, width="stretch", height=350)

# ---------------------------------------------------------------------------
# Tab: FFL results + null model
# ---------------------------------------------------------------------------

with tab_ffl:
    st.subheader("Feedforward loops (FFLs)")
    with st.expander("ℹ️ What's a feedforward loop, and what's a null model?", expanded=True):
        st.markdown(
            """
An **FFL** is a specific 3-gene wiring pattern: the master regulator WCC regulates
some other regulator *B*, and **both** WCC and *B* regulate the same target gene
*C*. It's counted here via the matrix formula the paper itself uses (the Hadamard
product `(R @ T) * T`, Section V-F) — this part is a direct reproduction of the
paper's method, not new.

What's new: the paper claims 71 FFLs is meaningful because it's more than a
textbook number from *E. coli* (40). That's a comparison to a *different network*
entirely. This tab instead asks: compared to a **randomized version of this exact
network** (same regulator target counts, same gene regulator-counts, edges
otherwise shuffled), is the observed FFL count unusually high? That's what the
histogram below shows.
"""
        )
    st.markdown(
        "An FFL is: WCC regulates a regulator *B*, and both WCC and *B* regulate the "
        "same target gene *C*. Counted via the Hadamard product `(R @ T) * T` "
        "(paper, Section V-F)."
    )
    col1, col2 = st.columns(2)
    col1.metric("Observed FFL count (this export)", ffl_result.total_ffl_count)
    col2.metric("Paper's reported FFL count", PAPER_REPORTED_FFL_COUNT)

    st.markdown("**FFL contribution per regulator (as the 'B' gene):**")
    contrib = ffl_result.per_regulator_ffl_count.sort_values(ascending=False)
    contrib.index = [f"{REGULATOR_INFO[lab]['gene']}" for lab in contrib.index]
    st.bar_chart(contrib)

    st.divider()
    st.subheader("Is the FFL count more than chance? (novel extension)")
    st.markdown(
        "A degree-preserving double-edge-swap null model: 1,000 randomized networks, "
        "each holding every regulator's target count and every gene's regulator-count "
        "fixed, re-counting FFLs each time."
    )

    if null_result is None:
        st.info(
            f"Null model not yet cached at `{NULL_CACHE_PATH}`. Run "
            "`python src/null_models.py` (~8-9 minutes) to generate it."
        )
    else:
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Observed", null_result.observed_ffl_count)
        col2.metric("Null mean ± SD", f"{null_result.mean_null:.1f} ± {null_result.std_null:.1f}")
        col3.metric("Z-score", f"{null_result.z_score:.2f}")
        col4.metric("Two-sided p-value", f"{null_result.p_value_two_sided:.3g}")
        st.plotly_chart(null_distribution_figure(null_result), width="stretch")
        st.markdown(
            "**Finding:** the observed FFL count is statistically indistinguishable "
            "from a randomly rewired network with the same degree sequence. See the "
            "Methodology tab for why this is a genuinely different conclusion from "
            "the paper's own qualitative comparison."
        )

# ---------------------------------------------------------------------------
# Tab: enrichment tests
# ---------------------------------------------------------------------------

with tab_enrich:
    st.subheader("Pairwise & triple-wise regulator co-targeting enrichment")
    with st.expander("ℹ️ What does this test mean, and how do I use the filters?", expanded=True):
        st.markdown(
            """
For every pair (55) and every group of three (165) of the 11 regulators, this asks:
*do they share more — or fewer — target genes than you'd expect by chance, given
how many targets each one has?* The chance baseline is an exact **hypergeometric
test** (the same math behind "what are the odds of drawing this many red balls from
an urn").

Because 220 tests are run at once, a raw p-value of 0.05 doesn't mean what it
usually means — by chance alone, ~11 of 220 tests would look "significant." The
**q-value** column fixes this (Benjamini-Hochberg FDR correction): filter on
`q-value < 0.05` and at most 5% of *those* results are expected to be false
positives, even across all 220 tests together.

**Using the filters below:** narrow to just pairs or just triples, just enriched or
just depleted relationships, and toggle "significant only" to hide anything that
doesn't clear the q < 0.05 bar. Click a column header in the table to sort by it.
"""
        )
    st.markdown(
        "All C(11,2) = 55 regulator pairs and C(11,3) = 165 regulator triples, tested "
        "against an exact hypergeometric null (pairs) / exact convolved-hypergeometric "
        "null (triples), Benjamini-Hochberg FDR corrected across all 220 tests at "
        "α = 0.05. The paper never ran this test."
    )

    if enrichment_df is None:
        st.info(
            f"Enrichment results not yet cached at `{ENRICHMENT_CACHE_PATH}`. Run "
            "`python src/enrichment.py` (~1 minute) to generate it."
        )
    else:
        n_sig = int(enrichment_df["significant_bh"].sum())
        col1, col2, col3 = st.columns(3)
        col1.metric("Total tests", len(enrichment_df))
        col2.metric("Significant (BH-FDR < 0.05)", n_sig)
        col3.metric("...all enriched or depleted?", enrichment_df.loc[enrichment_df["significant_bh"], "direction"].nunique())

        col_a, col_b, col_c = st.columns(3)
        test_type_filter = col_a.multiselect(
            "Test type", ["pair", "triple"], default=["pair", "triple"],
            help="Pair = 2 regulators compared; triple = 3 regulators compared at once.",
        )
        direction_filter = col_b.multiselect(
            "Direction", ["enriched", "depleted"], default=["enriched", "depleted"],
            help="Enriched = share more targets than chance; depleted = share fewer.",
        )
        sig_only = col_c.checkbox(
            "Significant only (BH-FDR < 0.05)", value=True,
            help="Hide any row whose corrected q-value is 0.05 or higher.",
        )

        filtered = enrichment_df[
            enrichment_df["test_type"].isin(test_type_filter) & enrichment_df["direction"].isin(direction_filter)
        ]
        if sig_only:
            filtered = filtered[filtered["significant_bh"]]
        filtered = filtered.sort_values("q_value")

        st.dataframe(
            filtered[
                ["regulators", "test_type", "observed_overlap", "expected_overlap", "direction", "p_value", "q_value", "significant_bh"]
            ],
            width="stretch",
            height=450,
        )
        st.markdown(
            "**Finding:** every significant result is a *depletion* — no regulator "
            "pair or triple shares more targets than chance. See the Methodology tab "
            "for why."
        )

# ---------------------------------------------------------------------------
# Tab: network visualization (interactive)
# ---------------------------------------------------------------------------

with tab_viz:
    st.subheader("Hub structure — click through the network")
    with st.expander("ℹ️ How to use this diagram", expanded=True):
        st.markdown(
            """
This reproduces the paper's Figure 1A layout: **WCC at the center**, broadcasting
out to the 11 downstream regulators (5 green transcriptional regulators, 6 orange
RNA operons). Node **size** = how many target genes that regulator has (its hub
size); the pale gray spokes are the paper's fixed WCC → regulator structure.

The colored **chords** between the outer regulators are this project's addition —
statistically significant co-targeting from the Enrichment Tests tab. Blue would
mean "these two regulators share more targets than chance" (enriched); red means
"fewer than chance" (depleted). *In this dataset, every significant relationship
turns out to be red/depleted — see the Methodology tab for why.*

**To explore:** pick a regulator from **"Focus on a regulator"** below to highlight
it and fade everything else into the background — its info, paper-sourced
description, and full list of significant relationships will appear beneath the
diagram. Use the **significance threshold** slider to loosen or tighten which
chords count as "significant," and the checkboxes to show only enriched or only
depleted relationships.
"""
        )

    focus_options = ["Show all regulators"] + net.regulator_names
    col_focus, col_thresh = st.columns([2, 1])
    focus_choice = col_focus.selectbox(
        "Focus on a regulator",
        options=focus_options,
        index=0,
        help="Highlights this regulator and its significant relationships; fades the rest.",
    )
    q_threshold = col_thresh.select_slider(
        "Significance threshold (q-value)",
        options=[0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0],
        value=0.05,
        help="A chord is drawn if its BH-adjusted q-value is at or below this. Loosen it to see more relationships.",
    )

    col_e, col_d = st.columns(2)
    show_enriched = col_e.checkbox("Show enriched (blue)", value=True, help="Regulator pairs sharing MORE targets than chance.")
    show_depleted = col_d.checkbox("Show depleted (red)", value=True, help="Regulator pairs sharing FEWER targets than chance.")

    focus_label = None
    if focus_choice != "Show all regulators":
        focus_idx = net.regulator_names.index(focus_choice)
        focus_label = net.regulator_labels[focus_idx]

    st.plotly_chart(
        hub_network_figure(
            net, enrichment_df,
            q_threshold=q_threshold, show_enriched=show_enriched, show_depleted=show_depleted,
            focus_label=focus_label,
        ),
        width="stretch",
    )

    if focus_label is not None:
        idx = net.regulator_labels.index(focus_label)
        with st.container(border=True):
            st.markdown(f"### {net.regulator_names[idx]} ({focus_label})")
            st.markdown(REGULATOR_BLURBS.get(focus_label, "_No description available._"))
            col1, col2, col3 = st.columns(3)
            col1.metric("Target genes (hub size)", int(net.T.loc[focus_label].sum()))
            col2.metric("Class", net.regulator_class[idx].replace("-", " "))
            col3.metric(
                "FFLs as the 'B' regulator",
                int(ffl_result.per_regulator_ffl_count.get(focus_label, 0))
                if focus_label != "Wcc" else "n/a (WCC is always 'A')",
            )

            if enrichment_df is not None:
                related = enrichment_df[
                    ((enrichment_df["reg_a"] == focus_label)
                     | (enrichment_df["reg_b"] == focus_label)
                     | (enrichment_df["reg_c"] == focus_label))
                    & (enrichment_df["q_value"] <= q_threshold)
                    & (enrichment_df["direction"].isin(
                        {d for d, show in [("enriched", show_enriched), ("depleted", show_depleted)] if show}
                    ))
                ].sort_values("q_value")
                st.markdown(f"**Significant relationships involving {net.regulator_names[idx]} at q ≤ {q_threshold}:**")
                if len(related):
                    st.dataframe(
                        related[["regulators", "test_type", "observed_overlap", "expected_overlap", "direction", "q_value"]],
                        width="stretch",
                        height=min(300, 45 + 35 * len(related)),
                    )
                else:
                    st.caption("None at this threshold — try loosening the significance slider above.")

    st.divider()
    st.subheader("Degree distribution")
    st.caption(
        "How many connections does a typical node have? Both regulators (many "
        "connections) and target genes (usually 1-2) are counted together, as in "
        "the paper's Figure 11."
    )
    st.markdown(
        f"Paper reports a power-law exponent of {PAPER_POWER_LAW_EXPONENT} across all "
        "network nodes (11 regulators + target genes)."
    )
    st.plotly_chart(degree_distribution_figure(stats_result), width="stretch")

    col1, col2, col3 = st.columns(3)
    col1.metric("Fitted exponent", f"{stats_result.power_law_exponent:.2f}", f"paper: {PAPER_POWER_LAW_EXPONENT}")
    col2.metric("Avg. node degree", f"{stats_result.average_node_degree:.2f}", f"paper: {PAPER_AVG_NODE_DEGREE}")
    col3.metric(
        "Singly / multiply regulated",
        f"{stats_result.singly_regulated} / {stats_result.multiply_regulated}",
        f"paper: {PAPER_SINGLY_REGULATED} / {PAPER_MULTIPLY_REGULATED}",
    )

# ---------------------------------------------------------------------------
# Tab: gene families (sequence-based, novel extension using the FASTA)
# ---------------------------------------------------------------------------

with tab_families:
    st.subheader("Sequence-based paralog gene families")
    with st.expander("ℹ️ What is this tab, and where does the data come from?", expanded=True):
        st.markdown(
            """
Everything else in this app comes from the paper's fitted network. This tab is
built from a completely different source: `data/raw/uniprotkb_proteome_UP000001805_*.fasta`,
the full UniProt reference proteome for *N. crassa* (10,255 protein sequences). The
paper never touches gene sequence data at all — this tab asks a question it can't:
**which clock-network genes have close evolutionary relatives (paralogs) elsewhere in
the genome, and how large are those paralog families?**

This matters for a companion project (the Scruse/Arnold/Robinson Partial-Duplication
framework), which needs real gene-duplication families to define an inheritance
probability. That framework couldn't be applied to this network before, for lack of
real sequence data — this FASTA is what unlocks it.

**Method, in two stages** (no BLAST/DIAMOND/MMseqs available in this environment):
1. **Cheap prefilter** — every protein is broken into overlapping 5-amino-acid
   fragments ("5-mers"); pairs of proteins sharing enough 5-mers become *candidates*.
   This is just a fast way to avoid comparing all 10,255 × 10,255 possible pairs.
2. **Real confirmation** — every candidate pair is then actually sequence-aligned
   (Biopython, BLOSUM62 substitution matrix, the same scoring approach BLAST itself
   uses) to get a real percent-identity. A pair only counts as "the same family" if
   that real alignment clears an identity and coverage bar — not the cheap prefilter.

A **family** is a group of genes connected (directly or transitively) by one of these
confirmed alignments. Family **size (k)** is just how many genes are in that group.
This tab specifically surfaces families of size **k=3** and **k=4**, plus the fuller
picture of every family that includes a clock-network gene.
"""
        )

    if families_df is None:
        st.info(
            f"Gene families not yet cached at `{FAMILIES_CACHE_PATH}`. Run "
            "`python src/sequence_families.py` (~3 minutes) to generate it."
        )
    else:
        st.markdown("### Regulator identifier resolution")
        st.caption(
            "Before building families, each of the paper's 11 regulators had to be "
            "matched to a protein sequence in the FASTA. This wasn't always clean — "
            "shown here exactly as found, not smoothed over."
        )
        res_rows = []
        for label, info in REGULATOR_INFO.items():
            ncu = info["ncu"]
            paper_name = info["gene"]
            resolved_name = regulator_gene_names[label]
            if resolved_name is None:
                reason = UNRESOLVED_REGULATORS.get(ncu) or EXCLUDED_REGULATORS.get(label, "unresolved")
                status = "Excluded" if label in EXCLUDED_REGULATORS else "Unresolved"
                res_rows.append((paper_name, ncu or "n/a", status, reason))
            else:
                via = f"alias -> {resolved_name}" if ncu in REGULATOR_ALIASES else "direct match"
                fam_row = families_df[families_df["gene_name"] == resolved_name]
                fam_size = int(fam_row["family_size"].iloc[0]) if len(fam_row) else 1
                res_rows.append((paper_name, ncu or "n/a", f"Resolved ({via})", f"family size {fam_size}"))
        res_df = pd.DataFrame(res_rows, columns=["Regulator", "NCU ID", "Status", "Detail"])
        st.dataframe(res_df, width="stretch", hide_index=True)

        n_resolved = sum(1 for r in res_rows if r[2].startswith("Resolved"))
        st.caption(
            f"{n_resolved} of 11 regulators resolved to a sequence. WCC is excluded "
            "(it's the WC-1/WC-2 protein complex, not a single gene product). "
            "NCU07392 (\"pro-1\" per the paper) is excluded because UniProt's own "
            "GN=pro-1 entry is a well-established, unrelated proline-biosynthesis "
            "enzyme — likely a gene-symbol collision, not the same locus; flagged "
            "rather than guessed. NCU03363 has no matching entry in this proteome at "
            "all, possibly a locus renumbered since the paper's 2022 gene models."
        )

        st.divider()
        st.markdown("### Family size distribution")
        size_counts = families_df.drop_duplicates("family_id")["family_size"].value_counts().sort_index()
        st.bar_chart(size_counts)
        st.caption(
            f"{int(size_counts.get(1, 0))} singleton genes with no "
            "detected paralog in this analysis; the rest fall into families of 2 or more."
        )

        st.divider()
        st.markdown("### Families of size k=3 and k=4")
        st.caption(
            "The specific question this tab was built to answer. Rows are grouped by "
            "family; genes already part of the clock network (paper's regulators or "
            "target genes) are marked."
        )
        k34 = families_df[families_df["family_size"].isin([3, 4])].sort_values(
            ["family_size", "family_id", "gene_name"]
        )
        if len(k34):
            for fam_id, group in k34.groupby("family_id"):
                k = int(group["family_size"].iloc[0])
                touches_clock = group["is_clock_network_gene"].any()
                badge = " 🔗 includes a clock-network gene" if touches_clock else ""
                with st.container(border=True):
                    st.markdown(f"**Family (k={k}){badge}**")
                    st.dataframe(
                        group[["gene_name", "description", "is_clock_network_gene"]],
                        width="stretch",
                        hide_index=True,
                    )
        else:
            st.caption("No k=3 or k=4 families found at current thresholds.")

        st.divider()
        st.markdown("### All clock-network genes with a detected paralog (any family size)")
        clock_multi = families_df[
            (families_df["is_clock_network_gene"]) & (families_df["family_size"] > 1)
        ].sort_values(["family_size", "family_id"], ascending=[False, True])
        st.caption(
            f"{clock_multi['family_id'].nunique()} families of size > 1 include at least one "
            "clock-network gene (regulator or target)."
        )
        st.dataframe(
            clock_multi[["gene_name", "family_id", "family_size", "description"]],
            width="stretch",
            height=300,
        )

        st.divider()
        st.markdown("### Look up any gene's family")
        lookup_options = sorted(families_df["gene_name"].unique())
        lookup_gene = st.selectbox(
            "Gene name or NCU ID",
            options=lookup_options,
            index=None,
            placeholder="Start typing…",
            help="Search across all 9,747 canonical genes in this analysis, not just clock-network ones.",
        )
        if lookup_gene:
            fam_id = families_df.loc[families_df["gene_name"] == lookup_gene, "family_id"].iloc[0]
            members = families_df[families_df["family_id"] == fam_id]
            st.dataframe(
                members[["gene_name", "description", "is_clock_network_gene"]],
                width="stretch",
                hide_index=True,
            )

        st.divider()
        with st.expander("⚠️ Limitations of this method"):
            st.markdown(
                """
- **k-mer prefilter can miss remote homologs.** Two genuinely related proteins that
  have diverged a lot can share too few exact 5-letter fragments to even become a
  candidate pair — this method has no gapped-alignment or profile/HMM step at the
  prefilter stage, so distant relationships can be invisible to it entirely.
- **Single-linkage clustering can over-merge.** If gene A resembles gene B, and gene B
  (separately) resembles unrelated gene C through a different shared feature (e.g. a
  common domain), all three can end up lumped into one family even if A and C aren't
  really related. Real orthology/family tools (e.g. OrthoFinder, MCL clustering) guard
  against this more carefully than a plain union-find does.
- **No formal statistical significance (no e-values).** Unlike BLAST, this doesn't
  compute the probability that a match this good could occur by chance — the identity/
  coverage thresholds are reasonable heuristics, not calibrated significance tests.
- **Same-gene duplicate UniProt records were explicitly filtered out** (pairs at ≥97%
  identity and ≥85% coverage are treated as the same physical gene recorded twice, not
  two paralogs) — this was necessary because the raw proteome file turned out to
  contain the identical protein under two different gene symbols in multiple cases
  (e.g. *mtA-1* / *matA-1*, *npl4* / NCU02680).

Treat every family reported here as a first-pass, biologically plausible hypothesis —
not a definitive, publication-grade orthology call.
"""
            )

# ---------------------------------------------------------------------------
# Tab: methodology & paper comparison
# ---------------------------------------------------------------------------

with tab_method:
    st.subheader("Reproduced from the paper vs. novel in this project")
    st.markdown(
        """
**Reproduced from the paper:**
- R/T matrix construction (Section III: "the regulatory CCG proteins do not regulate
  each other" — R is trivial, WCC activates all 11 regulators, they don't regulate
  each other).
- The FFL-counting formula via the Hadamard product `(R @ T) * T` (Section V-F).
- The per-regulator "B-gene" FFL breakdown via inner products of the WCC row of T
  against each other row (explicitly stated in Section V-F).
- The degree-distribution / hub-size / singly-vs-multiply-regulated comparisons
  (Section V-E).

**Novel in this project (the paper never computed these):**
- A degree-preserving double-edge-swap null model giving an empirical Z-score and
  p-value for the FFL count (`src/null_models.py`).
- Exact hypergeometric pairwise and convolved-hypergeometric triple-wise regulator
  co-targeting enrichment tests across all 220 combinations, Benjamini-Hochberg
  FDR corrected at α = 0.05 (`src/enrichment.py`).
"""
    )

    st.subheader("Validation table: this export vs. the published paper")
    validation_rows = [
        ("rok-1 (NCU00919) target count", 54, int(net.T.loc["PostReg9 NCU00919"].sum()), "Exact match"),
        ("lhp-1 (NCU08295) target count", 768, int(net.T.loc["PostReg5 NCU08295"].sum()), "Close, not exact"),
        ("Singly-regulated genes", PAPER_SINGLY_REGULATED, stats_result.singly_regulated, "Close, not exact"),
        ("Multiply-regulated genes", PAPER_MULTIPLY_REGULATED, stats_result.multiply_regulated, "Diverges"),
        ("Total clock genes", PAPER_TOTAL_GENES, net.T.shape[1], "Close, not exact"),
        ("Max regulators on one gene", "up to 11 (implied)", int(stats_result.gene_degrees.max()), "Diverges"),
        ("Total FFLs", PAPER_REPORTED_FFL_COUNT, ffl_result.total_ffl_count, "Diverges"),
        ("Genes with overlapping FFLs", 5, len(ffl_result.overlap_genes), "Diverges"),
        ("Avg. node degree", PAPER_AVG_NODE_DEGREE, round(stats_result.average_node_degree, 2), "Close"),
        ("Power-law exponent", PAPER_POWER_LAW_EXPONENT, round(stats_result.power_law_exponent, 2), "Close"),
    ]
    validation_df = pd.DataFrame(validation_rows, columns=["Quantity", "Paper", "This export", "Match?"])
    validation_df["Paper"] = validation_df["Paper"].astype(str)
    validation_df["This export"] = validation_df["This export"].astype(str)
    st.table(validation_df)

    st.subheader("Root cause of the divergence")
    st.markdown(
        """
Every gene in this exported matrix is regulated by **at most 2** of the 11
regulators. The paper's model assigns each target gene a continuous *binding
strength* per regulator (its Table 1; strengths across regulators sum to 1 per
gene, sampled via MCMC over the full ensemble). `Connection_Matrix.xlsx` looks
like a thresholded/binarized snapshot of that ensemble (e.g. top-1 or top-2
dominant regulator per gene) rather than the full multi-way soft assignment. That
single structural fact explains essentially every downstream divergence above: a
lower total FFL count, zero genes with overlapping FFLs, and — for the enrichment
tests — every significant result being a *depletion* (any two regulators' target
sets are pushed toward mutual exclusivity purely as a side effect of the
binarization, not necessarily real biological antagonism).

This is flagged rather than silently reconciled: report discrepancies, don't force
numbers to match.
"""
    )

    st.subheader("Citation")
    st.markdown(PAPER_CITATION)
