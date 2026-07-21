"""Interactive dashboard for the Neurospora crassa clock network project.

Reads the fitted regulator -> target matrix from Al-Omari, Griffith, Scruse,
Robinson, Schuttler & Arnold (2022), IEEE Access 10:32510-32524, and adds
statistical tests (a degree-preserving null model for the FFL count;
BH-FDR-corrected pairwise/triple/quartet regulator co-targeting enrichment) that the
paper itself did not perform.

Run with: streamlit run app/streamlit_app.py
"""
from __future__ import annotations

import base64
import os
import sys

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

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
from viz import (  # noqa: E402
    REGULATOR_COLORS,
    binding_strength_heatmap,
    degree_distribution_figure,
    ffl_network_figure,
    hub_network_figure,
    hub_network_figure_3d,
    null_distribution_figure,
    regulator_gene_figure_3d,
)
from link_analysis import build_elements, build_styles  # noqa: E402
from gene_layout_3d import (  # noqa: E402
    DEFAULT_CACHE_PATH as GENE_LAYOUT_3D_CACHE_PATH,
    load_result as load_gene_layout_3d,
)
from binding_strength import (  # noqa: E402
    DEFAULT_PATH as BINDING_STRENGTH_PATH,
    load_binding_strength_by_function,
)
from sequence_families import (  # noqa: E402
    DEFAULT_OUTPUT_PATH as FAMILIES_CACHE_PATH,
    EXCLUDED_REGULATORS,
    REGULATOR_ALIASES,
    UNRESOLVED_REGULATORS,
    resolve_regulator_gene_names,
)
from mutheta_loader import (  # noqa: E402
    DEFAULT_PATH as MUTHETA_PATH,
    N_MU,
    N_THETA,
    best_theta_to_regulator_mapping,
    compare_theta_dominance_to_real_network,
    load_mutheta_sim_with_genes,
)
from st_link_analysis import st_link_analysis  # noqa: E402

DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "Connection_Matrix.xlsx")

_DOCS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs")
FULL_PAPER_MD_PATH = os.path.join(_DOCS_DIR, "full_methodology_paper.md")
FULL_PAPER_PDF_PATH = os.path.join(_DOCS_DIR, "full_methodology_paper.pdf")
FULL_PAPER_DOCX_PATH = os.path.join(_DOCS_DIR, "full_methodology_paper.docx")

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

# One-sentence description of what each source file does, for anyone reading this
# app who doesn't have (or doesn't want to dig through) the underlying Python
# files. Used by file_ref() everywhere the app names a specific script.
FILE_BLURBS = {
    "data_loader.py": (
        "loads the raw Excel file, strips whitespace, resolves pandas' "
        "duplicate-column renaming, and separates the 3 core clock genes from "
        "the ~3,000 downstream target genes."
    ),
    "network_builder.py": (
        "builds the R (WCC → regulator) and T (regulator → target) matrices "
        "that every other script in this project starts from."
    ),
    "ffl_analysis.py": (
        "counts feedforward loops using the paper's own Hadamard-product "
        "formula, `(R @ T) * T`."
    ),
    "null_models.py": (
        "builds 1,000 randomly rewired versions of the network (each keeping "
        "every regulator's target count and every gene's regulator-count fixed) "
        "and recounts feedforward loops in each one, to test whether the real "
        "count is more than chance would predict."
    ),
    "enrichment.py": (
        "runs the hypergeometric (pairs) and convolved-hypergeometric (triples, "
        "quartets) significance tests for every regulator combination, then "
        "applies the Benjamini-Hochberg FDR correction across all 550 of them."
    ),
    "network_stats.py": (
        "computes the degree distribution, hub sizes, power-law fit, and the "
        "singly-/multiply-regulated gene split."
    ),
    "viz.py": (
        "builds the Plotly figures used throughout this app: the hub-and-spoke "
        "network diagram, the degree-distribution plot, and the null-model "
        "histogram."
    ),
    "link_analysis.py": (
        "builds the Cytoscape.js graph (elements + per-regulator color styles) "
        "for the 2D force-directed Regulators + Genes view."
    ),
    "gene_layout_3d.py": (
        "computes and caches a 3D force-directed layout (networkx "
        "`spring_layout(dim=3)`) for the same graph's 3D view."
    ),
    "binding_strength.py": (
        "loads the user-supplied binding-strength-by-function table and "
        "explains why its values are joint VTENS/MCMC point estimates, not "
        "independently testable per-parameter distributions."
    ),
    "sequence_families.py": (
        "parses the UniProt FASTA, finds candidate similar protein pairs via a "
        "cheap k-mer prefilter, confirms them with real BLOSUM62 sequence "
        "alignment, collapses duplicate UniProt records of the same gene, and "
        "groups the rest into paralog families."
    ),
}


def file_ref(filename: str) -> str:
    """Render a `src/filename.py` reference with its one-sentence description,
    so a reader without access to the Python source still knows what it does."""
    desc = FILE_BLURBS.get(filename)
    return f"`src/{filename}` (which {desc})" if desc else f"`src/{filename}`"


def tab_jump_button(
    label: str, key: str, target_tab_text: str = "Methodology", anchor_id: str | None = None
) -> None:
    """Render a button that switches to the top-level tab whose label contains
    `target_tab_text` (and, if given, scrolls to `anchor_id` within it).
    Streamlit's st.tabs has no built-in programmatic-switch API, so this
    clicks the tab's own button via a small injected script -- the only
    reliable way to jump across tabs in Streamlit.
    """
    if st.button(label, key=key):
        scroll_js = (
            f"""
            setTimeout(() => {{
                const el = doc.getElementById("{anchor_id}");
                if (el) el.scrollIntoView({{behavior: "smooth", block: "start"}});
            }}, 350);
            """
            if anchor_id
            else ""
        )
        components.html(
            f"""
            <script>
                const doc = window.parent.document;
                const tabs = doc.querySelectorAll('button[role="tab"]');
                for (const t of tabs) {{
                    if (t.innerText && t.innerText.includes("{target_tab_text}")) {{
                        t.click();
                        break;
                    }}
                }}
                {scroll_js}
            </script>
            """,
            height=0,
        )


# Backwards-compatible alias -- existing call sites jump to the Methodology tab.
def methodology_jump_button(label: str, key: str, anchor_id: str | None = None) -> None:
    tab_jump_button(label, key, target_tab_text="Methodology", anchor_id=anchor_id)


st.set_page_config(page_title="Neurospora Clock Network", layout="wide", page_icon="🕐")

# Shared column_config dicts so every results table shows human-readable
# headers and sensibly formatted numbers, not raw snake_case column names.
# Passing a key that isn't in the displayed DataFrame is harmless -- Streamlit
# just ignores it -- so the same dict is reused across tables with different
# column subsets.
ENRICHMENT_COLUMN_CONFIG = {
    "regulators": "Regulators",
    "test_type": "Test type",
    "observed_overlap": st.column_config.NumberColumn("Observed overlap", format="%d"),
    "expected_overlap": st.column_config.NumberColumn("Expected overlap", format="%.1f"),
    "direction": "Direction",
    "p_value": st.column_config.NumberColumn("p-value", format="scientific"),
    "q_value": st.column_config.NumberColumn("q-value", format="scientific"),
    "significant_bh": st.column_config.CheckboxColumn("Significant (BH-FDR)"),
}

FAMILY_COLUMN_CONFIG = {
    "gene_name": "Gene",
    "description": "Description",
    "is_clock_network_gene": st.column_config.CheckboxColumn("Clock-network gene?"),
    "family_id": st.column_config.NumberColumn("Family ID", format="%d"),
    "family_size": st.column_config.NumberColumn("Family size (k)", format="%d"),
}

MUTHETA_COLUMN_CONFIG = {
    "gene": "Gene",
    "dominant_theta": "Dominant θ column",
    **{f"mu_{i}": st.column_config.NumberColumn(f"μ_{i}", format="%.2f") for i in range(N_MU)},
    **{f"theta_{i}": st.column_config.NumberColumn(f"θ_{i}", format="%.3f") for i in range(N_THETA)},
}


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


@st.cache_resource
def _load_link_analysis_elements(_net):
    return build_elements(_net)


@st.cache_resource
def _load_link_analysis_styles(_net):
    return build_styles(_net)


@st.cache_resource
def _load_gene_layout_3d():
    return load_gene_layout_3d()


@st.cache_data
def _load_binding_strength():
    if not os.path.exists(BINDING_STRENGTH_PATH):
        return None
    return load_binding_strength_by_function()


@st.cache_data
def _load_mutheta_sim():
    if not os.path.exists(MUTHETA_PATH):
        return None
    return load_mutheta_sim_with_genes()


@st.cache_data
def _load_mutheta_mapping_result():
    return best_theta_to_regulator_mapping(connection_matrix_path=DATA_PATH)


@st.cache_data
def _load_mutheta_comparison():
    return compare_theta_dominance_to_real_network(connection_matrix_path=DATA_PATH)


@st.cache_data
def _load_full_paper_summary():
    """Extracts just the "Executive Summary" section of the full paper's
    markdown source (between that heading and the next top-level heading) --
    the in-app tab shows this short summary, not the whole document; the
    whole document is what the PDF/Word downloads and "open in new tab" are
    for."""
    if not os.path.exists(FULL_PAPER_MD_PATH):
        return None
    with open(FULL_PAPER_MD_PATH, "r", encoding="utf-8") as f:
        text = f.read()
    marker = "# Executive Summary"
    start = text.find(marker)
    if start == -1:
        return None
    start += len(marker)
    end = text.find("\n# ", start)
    section = text[start:end] if end != -1 else text[start:]
    section = section.strip()
    if section.endswith("---"):
        section = section[:-3].strip()
    return section


@st.cache_data
def _load_file_bytes(path: str):
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return f.read()


net = _load_network()
cleaning_report = _load_cleaning_report()
enrichment_df = _load_enrichment()
null_result = _load_null_model()
ffl_result = _load_ffl_result(net)
stats_result = _load_network_stats(net)
families_df = _load_gene_families()
link_elements = _load_link_analysis_elements(net)
link_node_styles, link_edge_styles = _load_link_analysis_styles(net)
gene_layout_3d_result = _load_gene_layout_3d()
binding_strength_df = _load_binding_strength()
mutheta_df = _load_mutheta_sim()
regulator_gene_names = resolve_regulator_gene_names(REGULATOR_INFO)
full_paper_summary = _load_full_paper_summary()
full_paper_pdf_bytes = _load_file_bytes(FULL_PAPER_PDF_PATH)
full_paper_docx_bytes = _load_file_bytes(FULL_PAPER_DOCX_PATH)

# ---------------------------------------------------------------------------
# Header (kept short so the tabs are visible without scrolling)
# ---------------------------------------------------------------------------

st.title("🕐 The Neurospora crassa Clock Network")
st.caption(
    "A statistical extension of the fitted RNA-operon / transcriptional-regulon "
    "network from Al-Omari et al. (2022), *IEEE Access*. New here? Start with the "
    "**Overview** tab."
)

tab_overview, tab_matrix, tab_binding, tab_mutheta, tab_ffl, tab_enrich, tab_viz, tab_families, tab_method, tab_paper = st.tabs(
    [
        "🏠 Overview",
        "📋 Matrix Explorer",
        "🔬 Binding Strengths",
        "🎲 Mu/Theta Sim",
        "🔁 FFLs & Null Model",
        "🧪 Enrichment Tests",
        "🕸️ Network Viz",
        "🧬 Gene Families",
        "📖 Methodology & Paper Comparison",
        "📄 Full Paper",
    ]
)

# ---------------------------------------------------------------------------
# Tab: Overview
# ---------------------------------------------------------------------------

with tab_overview:
    st.markdown(
        f"""
*Neurospora crassa* (bread mold) is the model organism used to work out how the
circadian clock operates at the molecular level. Its clock is driven by a core
feedback loop (*frq* / *wc-1* / *wc-2*), whose protein product **WCC** (White
Collar Complex) is the master regulator broadcasting time-of-day and light
information outward to the rest of the genome. This app explores the genome-scale
network **WCC** sits atop — fitted by *{PAPER_CITATION}* — and adds statistical
tests the paper itself never ran.
"""
    )

    col1, col2, col3 = st.columns(3)
    col1.metric("Regulators downstream of WCC", "11", help="5 transcriptional + 6 post-transcriptional (\"RNA operon\") regulators.")
    col2.metric("Target genes (cleaned)", f"{cleaning_report.n_target_genes_final:,}")
    col3.metric("Largest single hub", "lhp-1", help="768 targets per the paper -- more than the master regulator WCC itself.")

    with st.expander("📖 What this project adds (novel statistical tests)"):
        st.markdown(
            """
The fitted network is treated here as real input data (`data/raw/Connection_Matrix.xlsx`).
The paper reported network-motif counts (e.g. "71 feedforward loops, versus 40 in
*E. coli*") as a **qualitative** size comparison. This project asks the question the
paper didn't: *compared to what, exactly?* It adds a formal **degree-preserving null
model** (is the feedforward-loop count more than you'd expect from a randomly
rewired network with the same regulator target counts?), **Benjamini-Hochberg
FDR-corrected hypergeometric tests** for every pairwise and triple-wise regulator
combination (do any two or three regulators share more — or fewer — target genes
than chance?), and two data sources the paper never touches at all: real protein
sequence data (paralog gene families) and a simulated Mu/Theta dataset.
"""
        )

    with st.expander("⚠️ One important caveat, upfront"):
        st.markdown(
            """
The exported connection matrix this app reads turns out to cap every target gene at
**at most 2** regulators, whereas the paper's underlying model (**Variable Topology
Ensemble Methods**, VTENS) fits a continuous binding strength for every regulator
against every candidate target, with topology and kinetics estimated **jointly** by
MCMC under a "supernet" formulation — each accepted sample is one whole candidate
network, not an independent draw per (regulator, target) pair.
`Connection_Matrix.xlsx` is a thresholded snapshot of that joint fit (the dominant
regulator(s) per gene), not the full multi-way estimate. That single structural
difference explains essentially every place the numbers in this app diverge from
the paper's published figures. Every such divergence is reported explicitly rather
than adjusted away — see the **🔬 Binding Strengths** tab for a real
(function-level) slice of the underlying continuous estimates, and the **📖
Methodology & Paper Comparison** tab for the full accounting.
"""
        )

    st.divider()
    st.markdown("**How to use this app:**")
    st.markdown(
        """
- **📋 Matrix Explorer** — browse the cleaned data itself, and look up any gene.
- **🔬 Binding Strengths** — the paper's continuous, function-level binding strengths.
- **🎲 Mu/Theta Sim** — a simulated dataset joined to this project's genes, and why
  its leading interpretation didn't hold up.
- **🔁 FFLs & Null Model** — a 3-gene network pattern, tested against random chance.
- **🧪 Enrichment Tests** — do regulator pairs/triples share more or fewer targets
  than chance? A filterable results table.
- **🕸️ Network Viz** — interactive 2D/3D diagrams of the whole network.
- **🧬 Gene Families** — clock-network genes' paralogs, from real sequence data.
- **📖 Methodology & Paper Comparison** — every number here vs. the paper's, and
  why they differ.
"""
    )

# ---------------------------------------------------------------------------
# Tab: Matrix explorer
# ---------------------------------------------------------------------------

with tab_matrix:
    st.subheader("Raw connection matrix, cleaned")
    with st.expander("ℹ️ How to read this tab"):
        st.markdown(
            """
This tab shows the data behind everything else in the app: which of the 11
regulators targets which of the ~3,000 clock-controlled genes.

- The **raw vs. cleaned** table below shows how many targets each regulator has,
  before and after fixing whitespace/duplicate-column issues in the original file
  (exactly what "cleaning" means is spelled out just below, in its own dropdown).
- Use the **gene lookup dropdown** to pick any target gene and see which
  regulator(s) control it.
- The **full matrix** at the bottom is the actual 11 x ~3,000 binary table (1 =
  "this regulator targets this gene") that every other tab is computed from.
"""
        )

    with st.expander("🔧 What exactly does \"cleaning\" mean here? (the exact rule applied)"):
        st.markdown(
            f"""
The raw Excel export (`data/raw/Connection_Matrix.xlsx`) is **{cleaning_report.raw_shape[0]}
regulators × {cleaning_report.n_raw_columns} columns**, and isn't usable as-is. The
cleaning rule applied — in plain terms, in the order it's actually applied — is:

1. **Strip stray whitespace** from row and column labels (one row label and one
   column label in the raw file had trailing spaces that would otherwise make them
   look like different genes than they really are).
2. **Un-mangle duplicate columns.** When a spreadsheet has the same column header
   appear more than once, the tool that reads it (pandas) silently renames the
   repeats by appending `.1`, `.2`, etc. — so what looks like different genes
   (e.g. `NCU08295` and `NCU08295.1`) are actually the *same* gene column, listed
   twice in the original file. This step reverses that renaming and treats them as
   one column again. It found **{cleaning_report.n_duplicated_base_names} genes**
   that were duplicated this way, spread across
   **{cleaning_report.n_columns_in_duplicate_groups} raw columns**
   ({cleaning_report.n_unique_base_names} distinct genes remained after undoing it).
3. **Collapse duplicates.** For any gene that appeared in more than one (now
   un-mangled) column, a regulator is counted as targeting that gene if it shows a
   "1" in *any* of the duplicate columns (a logical OR).
4. **Split out the 3 core clock genes** (*wc-1*, *wc-2*, *frq*) into their own group,
   separate from the ~3,000 downstream target genes. These aren't "targets" in the
   paper's model — they're the clock mechanism itself, taken as already known and
   given, not something the regulators discovered/switched on. This tab's target-gene
   counts and dropdown exclude these
   ({', '.join(cleaning_report.core_clock_base_names)}).

After all four steps: **{cleaning_report.n_target_genes_final} target genes** remain
in the cleaned matrix used everywhere else in this app.

The code that does exactly this (and nothing more) is `src/data_loader.py` — a script
whose only job is turning the raw spreadsheet into the clean binary regulator×gene
table every other file in this project starts from. It's read directly, not
summarized from memory: every number above comes from actually running that
cleaning step against the raw file.
"""
        )

    col1, col2, col3 = st.columns(3)
    col1.metric("Raw columns", cleaning_report.n_raw_columns)
    col2.metric("Cleaned target genes", cleaning_report.n_target_genes_final)
    col3.metric("Core clock genes excluded", cleaning_report.n_core_clock_columns)

    st.markdown("**Regulator target counts (row sums), raw vs. cleaned:**")
    cmp = pd.DataFrame({"raw": cleaning_report.raw_row_sums, "cleaned": cleaning_report.cleaned_row_sums})
    cmp.index = [f"{lab}  —  {REGULATOR_INFO.get(lab.strip(), {}).get('gene', '?')}" for lab in cmp.index]
    st.dataframe(
        cmp,
        width="stretch",
        column_config={"raw": "Raw target count", "cleaned": "Cleaned target count"},
    )

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
# Tab: binding strength by function (user-supplied continuous data)
# ---------------------------------------------------------------------------

with tab_binding:
    st.subheader("Binding strength by function")
    with st.expander("ℹ️ What is this, and what it is *not*"):
        st.markdown(
            f"""
A separate, user-supplied table (`data/raw/binding_strength_by_function.csv`):
21 KEGG/GO functional categories (ribosome, proteasome, carbon metabolism,
etc.) x the same 11 regulators, with a continuous value in each cell instead
of the binary 0/1 of the Matrix Explorer tab.

**What these numbers are not — read this before drawing any conclusion from
them.** It's tempting to read a `0.0000` cell as "this regulator's ensemble
distribution for this function excluded zero" and a positive cell as
"included it" — that would be a real statistical claim. **It is not what
this table supports.** The paper's Variable Topology Ensemble Method
(VTENS) does not fit each regulator-function binding strength as an
independent parameter with its own marginal posterior, later tested against
zero. Under VTENS's "supernet" formulation, every regulator is initially
allowed to regulate every candidate function, and the network *topology*
and its *kinetic parameters* are estimated **jointly** by MCMC — each
accepted sample is one entire candidate network consistent with the
underlying mRNA time-course data, not a draw from 231 independent
per-(function, regulator) distributions. The values below are that joint
ensemble's point estimates (the same estimates used to pick the dominant
regulator for each function/gene and build `Connection_Matrix.xlsx` in the
first place) — not the result of, and not re-derivable into, per-cell
hypothesis tests. A blank/zero cell means "not the estimated dominant
regulator for this function," full stop.
"""
        )

    if binding_strength_df is None:
        st.info(f"Binding-strength table not found at `{BINDING_STRENGTH_PATH}`.")
    else:
        st.plotly_chart(binding_strength_heatmap(net, binding_strength_df), width="stretch")

        st.divider()
        st.markdown("**Full table**")
        display_df = binding_strength_df.copy()
        display_df.columns = net.regulator_names
        st.dataframe(
            display_df,
            width="stretch",
            height=350,
            column_config={
                name: st.column_config.NumberColumn(name, format="%.4f")
                for name in display_df.columns
            },
        )

# ---------------------------------------------------------------------------
# Tab: simulated Mu/Theta dataset (user-supplied, gene identity confirmed,
# regulator-assignment hypothesis tested and refuted)
# ---------------------------------------------------------------------------

with tab_mutheta:
    st.subheader("Simulated Mu/Theta dataset")
    with st.expander("ℹ️ What is this, and what did we learn from it?"):
        st.markdown(
            """
Two more user-supplied files: `data/raw/MuThetaDataSim_m4.txt` (2,418 simulated
records, CR-delimited, no header) and `data/raw/ALLNCU.txt` (its gene IDs, supplied
afterward). Their row order is a confirmed one-to-one correspondence — this is the
indexing used throughout the analysis to associate each record's Mu/Theta values
with its gene, not just an inference from matching line counts.

Each record has two blocks:
- **Mu** (10 values) — continuous, roughly 0-100. **Confirmed:** these are
  *conditional* ensemble averages, not simple means, over an ensemble of 2,000
  MCMC runs. For each ensemble member, the regulator with the largest inferred
  binding strength for that gene is identified as dominant; a Mu column's value is
  the average binding strength of its regulator *only* over the ensemble members
  where that regulator was dominant for the gene — not an average over all 2,000
  members. Still open: which of this project's 10 non-WCC regulators each of
  `mu_0`..`mu_9` positionally corresponds to.
- **Theta** (6 values) — always sums to exactly 1.0: a probability/mixture-weight
  vector over 6 categories. The conditional-averaging clarification above is
  specific to Mu; Theta's identity is separately tested (and refuted, for the one
  hypothesis tried) below.

**Confirmed:** every one of the 2,216 distinct genes in `ALLNCU.txt` is a real
target gene in this project's `Connection_Matrix.xlsx` — this dataset shares its
gene space with the rest of the project, not a coincidence.

**Tested and refuted:** the leading hypothesis was that Theta's 6 categories are
each gene's soft assignment across this project's 6 real post-transcriptional
"RNA operon" regulators — motivated by the matching count of 6, and Theta's most
common category dominating similarly to how lhp-1 dominates real target counts.
Brute-forcing all 720 ways to label the 6 Theta columns against the 6 real
regulators, then scoring per-gene prediction accuracy (does argmax(Theta) match a
gene's real single post-transcriptional regulator?) against the 887 sim genes with
one unambiguous real regulator, shows the *best possible* labeling still loses to
simply guessing lhp-1 for every gene — see the metrics below. What Mu/Theta
actually represent remains open; full detail in the README.
"""
        )

    if mutheta_df is None:
        st.info(f"Mu/Theta data not found at `{MUTHETA_PATH}`.")
    else:
        mapping_result = _load_mutheta_mapping_result()
        majority_gene_name = REGULATOR_INFO[mapping_result["majority_label"]]["gene"]

        col1, col2, col3 = st.columns(3)
        col1.metric(
            "Sim records / distinct genes",
            f"{len(mutheta_df):,} / {mutheta_df['gene'].nunique():,}",
        )
        col2.metric(
            "Best-case accuracy (any labeling)",
            f"{mapping_result['best_acc']:.1%}",
            help="Best of all 720 permutations of theta column -> real regulator label.",
        )
        col3.metric(
            f"Majority baseline (always '{majority_gene_name}')",
            f"{mapping_result['majority_baseline_acc']:.1%}",
            delta=f"{(mapping_result['best_acc'] - mapping_result['majority_baseline_acc']):.1%} vs. best-case",
            delta_color="inverse",
            help="Always guessing the most common real regulator beats Theta's argmax under every labeling tried.",
        )
        st.caption(
            "Best-case accuracy is *below* the majority baseline — Theta's argmax does not "
            "predict a gene's real post-transcriptional regulator better than always "
            "guessing lhp-1. Hypothesis refuted, not confirmed."
        )

        st.divider()
        st.markdown("**Theta's dominant-category share: simulated vs. real regulator target share**")
        st.caption(
            "If Theta really encoded regulator assignment, these two bars would track each "
            "other closely for every regulator. Only the biggest one (lhp-1) roughly does."
        )
        comparison_df = _load_mutheta_comparison().copy()
        comparison_df["Regulator"] = [
            REGULATOR_INFO[lab]["gene"] for lab in comparison_df["candidate_regulator"]
        ]
        chart_df = comparison_df.set_index("Regulator")[["sim_dominant_share", "real_target_share"]]
        chart_df.columns = ["Simulated dominant share", "Real target share"]
        st.bar_chart(chart_df)

        st.divider()
        st.markdown("**Look up a gene's simulated Mu/Theta values**")
        mutheta_gene_options = sorted(mutheta_df["gene"].unique())
        mutheta_gene_query = st.selectbox(
            "Gene (NCU ID)",
            options=mutheta_gene_options,
            index=None,
            placeholder="Start typing an NCU ID…",
            help="Some genes appear more than once (202 repeats in ALLNCU.txt) — each has its own row.",
            key="mutheta_gene_lookup",
        )
        if mutheta_gene_query:
            rows = mutheta_df[mutheta_df["gene"] == mutheta_gene_query].reset_index()
            st.dataframe(
                rows,
                width="stretch",
                hide_index=True,
                column_config=MUTHETA_COLUMN_CONFIG,
            )

        st.divider()
        st.markdown("**Full simulated dataset**")
        st.caption("Rows = sim records, joined to gene names by row order. Scroll to explore.")
        st.dataframe(
            mutheta_df.reset_index(),
            width="stretch",
            height=350,
            hide_index=True,
            column_config=MUTHETA_COLUMN_CONFIG,
        )

# ---------------------------------------------------------------------------
# Tab: FFL results + null model
# ---------------------------------------------------------------------------

with tab_ffl:
    st.subheader("Feedforward loops (FFLs)")
    with st.expander("ℹ️ What's a feedforward loop, and what's a null model?"):
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
    st.caption(f"Counted here by {file_ref('ffl_analysis.py')}.")
    col1, col2 = st.columns(2)
    col1.metric("Observed FFL count (this export)", ffl_result.total_ffl_count)
    col2.metric("Paper's reported FFL count", PAPER_REPORTED_FFL_COUNT)

    st.markdown("**FFL contribution per regulator (as the 'B' gene):**")
    contrib = ffl_result.per_regulator_ffl_count.sort_values(ascending=False)
    contrib.index = [f"{REGULATOR_INFO[lab]['gene']}" for lab in contrib.index]
    contrib.index.name = "Regulator"
    contrib.name = "FFL count"
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
            f"`python src/null_models.py` (~8-9 minutes) to generate it — that's "
            f"{file_ref('null_models.py')}."
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
    st.subheader("Pairwise, triple-wise & quartet-wise regulator co-targeting enrichment")
    with st.expander("ℹ️ What does this test mean, and how do I use the filters?"):
        st.markdown(
            """
For every pair (55), every group of three (165), and every group of four (330) of
the 11 regulators, this asks: *do they share more — or fewer — target genes than
you'd expect by chance, given how many targets each one has?* The chance baseline
is an exact **hypergeometric test** (the same math behind "what are the odds of
drawing this many red balls from an urn"), extended to triples and quartets by
sequential convolution (see the module docstring in `enrichment.py`).

Because 550 tests are run at once, a raw p-value of 0.05 doesn't mean what it
usually means — by chance alone, ~28 of 550 tests would look "significant." The
**q-value** column fixes this (Benjamini-Hochberg FDR correction): filter on
`q-value ≤ 0.05` and at most 5% of *those* results are expected to be false
positives, even across all 550 tests together.

**Using the filters below:** narrow to pairs, triples, and/or quartets, just
enriched or just depleted relationships, and toggle "significant only" to hide
anything that doesn't clear the q ≤ 0.05 bar. Click a column header in the table
to sort by it.
"""
        )
    st.caption(
        "Pairs use an exact hypergeometric test; triples and quartets use a "
        "convolved-hypergeometric test (see module docstring for the derivation) — "
        f"all run by {file_ref('enrichment.py')}. The paper never ran any of them."
    )

    if enrichment_df is None:
        st.info(
            f"Enrichment results not yet cached at `{ENRICHMENT_CACHE_PATH}`. Run "
            f"`python src/enrichment.py` (~1 minute) to generate it — that's "
            f"{file_ref('enrichment.py')}."
        )
    else:
        n_sig = int(enrichment_df["significant_bh"].sum())
        sig_direction_counts = enrichment_df.loc[enrichment_df["significant_bh"], "direction"].value_counts()
        n_sig_enriched = int(sig_direction_counts.get("enriched", 0))
        n_sig_depleted = int(sig_direction_counts.get("depleted", 0))

        col1, col2, col3 = st.columns(3)
        col1.metric("Total tests", len(enrichment_df))
        col2.metric("Significant (BH-FDR ≤ 0.05)", n_sig)
        col3.metric(
            "...of which enriched vs. depleted",
            f"{n_sig_enriched} / {n_sig_depleted}",
            help=(
                "Of the significant results, how many are enrichments (more shared "
                "targets than chance) vs. depletions (fewer than chance). "
                f"{n_sig_enriched} enriched, {n_sig_depleted} depleted here — use "
                "the button below the table to jump to why depletions dominate so "
                "completely."
            ),
        )
        methodology_jump_button(
            "📖 Why depletions dominate — jump to Methodology",
            key="jump_methodology_top",
            anchor_id="why-depleted-anchor",
        )

        col_a, col_b, col_c = st.columns(3)
        test_type_filter = col_a.multiselect(
            "Test type", ["pair", "triple", "quartet"], default=["pair", "triple", "quartet"],
            help="Pair = 2 regulators compared; triple = 3; quartet = 4, all at once.",
        )
        direction_filter = col_b.multiselect(
            "Direction", ["enriched", "depleted"], default=["enriched", "depleted"],
            help="Enriched = share more targets than chance; depleted = share fewer.",
        )
        sig_only = col_c.checkbox(
            "Significant only (BH-FDR ≤ 0.05)", value=True,
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
            hide_index=True,
            column_config=ENRICHMENT_COLUMN_CONFIG,
        )
        st.markdown(
            "**Finding:** every significant result is a *depletion* — no regulator "
            "pair, triple, or quartet shares more targets than chance."
        )
        methodology_jump_button(
            "📖 See the Methodology tab for why", key="jump_methodology_bottom",
            anchor_id="why-depleted-anchor",
        )

# ---------------------------------------------------------------------------
# Tab: network visualization (interactive)
# ---------------------------------------------------------------------------

with tab_viz:
    subtab_hub, subtab_genes, subtab_ffl = st.tabs(
        ["🕸️ Hub structure", "🌐 Regulators + Genes", "🔁 FFL Network"]
    )

    with subtab_hub:
        st.subheader("Hub structure — click through the network")
        with st.expander("ℹ️ How to use this diagram"):
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

**2D vs. 3D:** both show the same data. In **3D**, the two regulator classes are
also separated vertically — the 5 green transcriptional regulators form a ring
above WCC, the 6 orange RNA operons a ring below it — so class membership reads
from height as well as color. Click and drag inside a 3D diagram to rotate it,
scroll (or pinch) to zoom.

**To explore:** pick a regulator from **"Focus on a regulator"** below to highlight
it and fade everything else into the background — its info, paper-sourced
description, and full list of significant relationships will appear beneath the
diagram. Use the **significance threshold** slider to loosen or tighten which
chords count as "significant," and the checkboxes to show only enriched or only
depleted relationships.
"""
            )
            methodology_jump_button(
                "📖 Why every chord is red — jump to Methodology",
                key="jump_methodology_viz",
                anchor_id="why-depleted-anchor",
            )

        view_mode = st.segmented_control(
            "Diagram view",
            options=["2D", "3D"],
            default="2D",
            help="3D separates the two regulator classes vertically and lets you drag to rotate.",
        )
        view_mode = view_mode or "2D"
        if view_mode == "3D":
            st.caption("Drag to rotate, scroll to zoom. Green ring above WCC = transcriptional; orange ring below = RNA operons.")

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

        figure_fn = hub_network_figure_3d if view_mode == "3D" else hub_network_figure
        st.plotly_chart(
            figure_fn(
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
                col2.metric("Class", net.regulator_class[idx].capitalize())
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
                            hide_index=True,
                            column_config=ENRICHMENT_COLUMN_CONFIG,
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

    with subtab_genes:
        st.subheader(f"All {net.T.shape[1]:,} target genes, hung off their regulators")

        genes_view_mode = st.segmented_control(
            "Diagram view",
            options=["2D", "3D"],
            default="2D",
            key="genes_view_mode",
            help="Same graph, same per-regulator colors, either engine: 2D uses Cytoscape.js's "
                 "fcose in the browser, 3D uses a precomputed networkx force-directed layout.",
        )
        genes_view_mode = genes_view_mode or "2D"

        if genes_view_mode == "2D":
            with st.expander("ℹ️ How to read this diagram"):
                st.markdown(
                    """
Every one of this project's 3,026 target genes plotted as its own point,
each connected by an edge to the regulator(s) that target it — laid out by
Cytoscape.js's **fcose** force-directed algorithm rather than a fixed ring,
so each regulator's targets pull together into their own starburst cluster
purely from the pull of their shared edges.

Each regulator gets **its own color** (11 total — see the legend below the
diagram), applied to that regulator's node, its exclusively-targeted genes,
and the edges between them, so a cluster's hub is identifiable on sight. The
genes targeted by **two** regulators (this export's maximum per gene) sit
between their two clusters in neutral gray, with an edge to each.

**Click a node** to highlight its direct neighbors and fade the rest —
built into the component, no extra controls needed. Drag to pan, scroll to
zoom, and use the toolbar in the top-right corner for fullscreen or a JSON
export of the graph.

**If the graph looks too zoomed-in when this tab first loads:** click the
target-shaped **"Fit to Screen"** icon in the top-right toolbar once — the
component doesn't always auto-fit its very first render to the full width.
"""
                )

            st_link_analysis(
                link_elements,
                # wider than fcose's defaults (nodeRepulsion=4500, idealEdgeLength=50) --
                # at the default spacing, regulator hubs and their gene clusters sit too
                # close together to see individual edges between them.
                layout={
                    "name": "fcose",
                    "padding": 20,
                    "animationDuration": 500,
                    "fit": True,
                    "animate": True,
                    "nodeDimensionsIncludeLabels": True,
                    "nodeRepulsion": 12000,
                    "idealEdgeLength": 120,
                },
                node_styles=link_node_styles,
                edge_styles=link_edge_styles,
                height=750,
                key="regulator_gene_link_analysis",
            )

            legend_cols = st.columns(6)
            for i, lab in enumerate(net.regulator_labels):
                with legend_cols[i % 6]:
                    color = REGULATOR_COLORS[lab]
                    st.markdown(
                        f"<span style='color:{color}'>&#9679;</span> {net.regulator_names[i]}",
                        unsafe_allow_html=True,
                    )
        else:
            with st.expander("ℹ️ How to read this diagram"):
                st.markdown(
                    """
The same regulator + gene graph as the 2D view, laid out in true 3D instead:
positions come from **networkx's `spring_layout(dim=3)`**, a real
force-directed physics simulation (not the fixed-ring layout the Hub
structure subtab uses) — the same idea as the 2D fcose graph, computed by a
different engine and precomputed/cached rather than run live in the browser
(it takes ~45s for this graph's 3,037 nodes).

Same per-regulator coloring as everywhere else in this tab. With ~3,026
gene points on screen, use **"Focus on a regulator"** below to isolate one
hub's genes at full opacity and fade the rest. Drag to rotate, scroll to
zoom.
"""
                )

            if gene_layout_3d_result is None:
                st.info(
                    f"3D layout not yet cached at `{GENE_LAYOUT_3D_CACHE_PATH}`. Run "
                    f"`python src/gene_layout_3d.py` (~45s) to generate it — that's "
                    f"{file_ref('gene_layout_3d.py')}."
                )
            else:
                focus_choice_genes_3d = st.selectbox(
                    "Focus on a regulator",
                    options=["Show all regulators"] + net.regulator_names,
                    index=0,
                    key="genes_3d_focus",
                    help="Isolate this regulator's target genes at full opacity; fade everything else.",
                )
                focus_label_genes_3d = None
                if focus_choice_genes_3d != "Show all regulators":
                    focus_label_genes_3d = net.regulator_labels[net.regulator_names.index(focus_choice_genes_3d)]

                st.plotly_chart(
                    regulator_gene_figure_3d(net, gene_layout_3d_result, focus_label=focus_label_genes_3d),
                    width="stretch",
                )

    with subtab_ffl:
        st.subheader(f"Feedforward-loop network ({ffl_result.total_ffl_count} loops)")
        with st.expander("ℹ️ How to read this diagram"):
            st.markdown(
                """
This is the same feedforward-loop (FFL) count from the **FFLs & Null Model**
tab — WCC (A) regulates a regulator (B), and both WCC and B regulate the same
target gene (C) — laid out here as the actual loops it's counting, rather than
just the summary numbers.

**WCC sits at the center**; every regulator that completes at least one FFL as
the "B" gene sits in a ring around it (gray spokes = the fixed WCC → B
structure the paper specifies). Each B's FFL-completing target genes hang off
it as small dots — the **solid colored line** is the B → C edge, the **dotted
violet line** is the WCC → C edge that closes the same triangle back at the
center.

Use **"Focus on a regulator"** to isolate one regulator's loops.
"""
            )

        ffl_focus_names = [n for n in net.regulator_names if n != "WCC"]
        focus_choice_ffl = st.selectbox(
            "Focus on a regulator",
            options=["Show all regulators"] + ffl_focus_names,
            index=0,
            key="ffl_focus",
            help="Isolate this regulator's feedforward loops; fade everything else.",
        )
        focus_label_ffl = None
        if focus_choice_ffl != "Show all regulators":
            focus_label_ffl = net.regulator_labels[net.regulator_names.index(focus_choice_ffl)]

        if ffl_result.total_ffl_count == 0:
            st.info("No feedforward loops in this export — nothing to draw.")
        else:
            st.plotly_chart(
                ffl_network_figure(net, ffl_result, focus_label=focus_label_ffl),
                width="stretch",
            )

# ---------------------------------------------------------------------------
# Tab: gene families (sequence-based, novel extension using the FASTA)
# ---------------------------------------------------------------------------

with tab_families:
    st.subheader("Sequence-based paralog gene families")
    with st.expander("ℹ️ What is this tab, and where does the data come from?"):
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
            f"`python src/sequence_families.py` (~3 minutes) to generate it — "
            f"that's {file_ref('sequence_families.py')}."
        )
    else:
        n_clock = int(families_df["is_clock_network_gene"].sum())
        n_total = len(families_df)
        n_not_clock = n_total - n_clock

        st.markdown("### What counts as a \"clock-network gene\" here — and why most genes don't")
        st.caption(
            f"Of the {n_total:,} genes in this analysis, only {n_clock:,} ({n_clock / n_total:.0%}) "
            f"are marked as clock-network genes. The other {n_not_clock:,} ({n_not_clock / n_total:.0%}) "
            "are real *N. crassa* genes too — just not part of the paper's clock model."
        )
        with st.expander("ℹ️ Why aren't most genes \"clock-network genes\"?"):
            st.markdown(
                f"""
**The definition used in this tab:** a gene counts as a "clock-network gene" if it's
either (a) one of the paper's 11 regulators, or (b) one of the
{cleaning_report.n_target_genes_final:,} target genes in the paper's fitted network
(`data/raw/Connection_Matrix.xlsx`, the same matrix explored in the Matrix Explorer
tab). Everything else in the FASTA — the other {n_not_clock:,} genes — is marked
`is_clock_network_gene = False`.

**Why the majority aren't clock-network genes, in order of how much each explains:**

1. **The FASTA is the *entire* genome; the clock network is a deliberately chosen
   slice of it.** `data/raw/uniprotkb_proteome_UP000001805_*.fasta` contains every
   protein-coding gene *N. crassa* has — roughly 10,000+ genes total, covering every
   biological process the organism does (metabolism, structure, reproduction, stress
   response, etc.), not just the clock.
2. **The paper itself only ever selected ~31% of the genome as "clock-controlled."**
   Per the paper's own numbers, of about 11,000 total *N. crassa* genes, only 3,380
   showed a circadian rhythm and/or a light-entrainment response in their RNA
   profiling experiments and could be confidently assigned to a regulator in their
   model fit. The other ~69% of the genome was never a candidate for "clock-network
   gene" status in the first place — those genes are presumably doing other jobs
   (housekeeping, metabolism, development, stress response) with no measured
   circadian signature at all.
3. **This app's own sequence-matching step lost a few more.** Of the
   {cleaning_report.n_target_genes_final:,} target genes and 11 regulators the paper
   *did* identify, only {n_clock:,} could actually be matched by name to a sequence in
   this specific UniProt download (see the regulator resolution table below, and the
   Matrix Explorer tab's cleaning-rule dropdown for the target-gene side of this same
   problem) — some genes are recorded under a different symbol than the one the paper
   used, and a few weren't found at all. So the count here is a slight *undercount* of
   the paper's true clock network, not just "proteome minus clock network."

In short: this isn't a data quality problem to fix — it's exactly what you'd expect,
since a circadian clock is one biological program among many, and most of an
organism's genome is doing something else entirely.
"""
            )

        st.divider()
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
        n_singletons = int(size_counts.get(1, 0))
        st.caption(
            f"{n_singletons:,} genes ({n_singletons / size_counts.sum():.0%} of all "
            f"{int(size_counts.sum()):,} families) have no detected paralog at all — "
            "so large a majority that including them here would make every other bar "
            "invisible. The chart below is families of size ≥ 2 only."
        )
        multi_size_counts = size_counts[size_counts.index > 1]
        multi_size_counts.index.name = "Family size (k)"
        multi_size_counts.name = "Number of families"
        st.bar_chart(multi_size_counts)

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
                        column_config=FAMILY_COLUMN_CONFIG,
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
            hide_index=True,
            column_config=FAMILY_COLUMN_CONFIG,
        )

        st.divider()
        st.markdown("### Look up any gene's family")
        lookup_options = sorted(families_df["gene_name"].unique())
        lookup_gene = st.selectbox(
            "Gene name or NCU ID",
            options=lookup_options,
            index=None,
            placeholder="Start typing…",
            help=(
                f"Search across all {len(lookup_options):,} canonical genes in this "
                "analysis, not just clock-network ones."
            ),
        )
        if lookup_gene:
            fam_id = families_df.loc[families_df["gene_name"] == lookup_gene, "family_id"].iloc[0]
            members = families_df[families_df["family_id"] == fam_id]
            st.dataframe(
                members[["gene_name", "description", "is_clock_network_gene"]],
                width="stretch",
                hide_index=True,
                column_config=FAMILY_COLUMN_CONFIG,
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
    with st.expander("📋 Full breakdown, file by file"):
        st.markdown(
            f"""
**Reproduced from the paper** — the code lives in {file_ref('network_builder.py')}
and {file_ref('ffl_analysis.py')}:
- R/T matrix construction (Section III: "the regulatory CCG proteins do not regulate
  each other" — R is trivial, WCC activates all 11 regulators, they don't regulate
  each other).
- The FFL-counting formula via the Hadamard product `(R @ T) * T` (Section V-F).
- The per-regulator "B-gene" FFL breakdown via inner products of the WCC row of T
  against each other row (explicitly stated in Section V-F).
- The degree-distribution / hub-size / singly-vs-multiply-regulated comparisons
  (Section V-E), computed by {file_ref('network_stats.py')}.
- The function-level continuous binding-strength table (Binding Strengths tab) is
  the paper's own joint VTENS/MCMC point estimates, transcribed as supplied —
  {file_ref('binding_strength.py')}.

**Novel in this project (the paper never computed these):**
- A degree-preserving double-edge-swap null model giving an empirical Z-score and
  p-value for the FFL count — {file_ref('null_models.py')}.
- Exact hypergeometric pairwise and convolved-hypergeometric triple- and
  quartet-wise regulator co-targeting enrichment tests across all 550 combinations,
  Benjamini-Hochberg FDR corrected at α = 0.05 — {file_ref('enrichment.py')}.
- Sequence-based paralog gene families from real protein sequence data —
  {file_ref('sequence_families.py')}.
- Force-directed layouts of every regulator plus all 3,026 target genes — the
  paper's own Figure 1A is a fixed hub-and-spoke — in 2D via Cytoscape.js
  ({file_ref('link_analysis.py')}) and true 3D via networkx
  ({file_ref('gene_layout_3d.py')}).

All charts in every tab of this app are built by {file_ref('viz.py')}, shared
across tabs so the same colors and layout logic are used everywhere.
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

    st.markdown('<div id="why-depleted-anchor"></div>', unsafe_allow_html=True)
    st.subheader("Root cause of the divergence")
    with st.expander("🔍 Why the numbers don't match", expanded=True):
        st.markdown(
            """
Every gene in this exported matrix is regulated by **at most 2** of the 11
regulators. The paper's underlying model (VTENS) doesn't fit each gene a set of
11 independent per-regulator binding strengths tested one at a time — under its
"supernet" formulation, every regulator is initially allowed to regulate every
candidate target, and the network *topology* and its *kinetic parameters* are
estimated **jointly** by MCMC: each accepted sample is one entire candidate
network consistent with the mRNA time-course data, not 11 independent draws per
gene. `Connection_Matrix.xlsx` is a thresholded snapshot of that joint fit (the
dominant regulator(s) per gene) rather than the full multi-way estimate — see the
**🔬 Binding Strengths** tab for a real (function-level) look at what the
continuous estimates look like and what they don't support. That single
structural fact explains essentially every downstream divergence above: a lower
total FFL count, zero genes with overlapping FFLs, and — for the enrichment tests
— every significant result being a *depletion* (any combination of regulators'
target sets is pushed toward mutual exclusivity purely as a side effect of the
binarization, not necessarily real biological antagonism).

This is flagged rather than silently reconciled: report discrepancies, don't force
numbers to match.
"""
        )

    st.subheader("Citation")
    st.markdown(PAPER_CITATION)

    st.divider()
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
  once (550, here), so "5% of my *significant* results are false positives" stays
  true even after testing every regulator pair, triple, and quartet. Always use
  q-value, not raw p-value, to decide significance when many tests are run together.
- **Enriched vs. depleted** — enriched = two regulators share *more* targets than
  chance predicts; depleted = *fewer* than chance predicts.
- **Paralog / gene family** — genes descended from a common ancestral gene via
  duplication, still recognizably similar in sequence. A gene "family" here is a
  group of such genes; family *size* (k) is how many genes are in it.
- **Canonical gene** — after collapsing duplicate UniProt records of the same
  physical gene (see Gene Families tab), the one name kept to represent it.
"""
        )

    st.divider()
    st.markdown(
        "**Want the full story, written for a general audience?** The **Full Paper** "
        "tab has a complete, plain-language write-up of everything in this app — the "
        "biology, every method, every result, and why the numbers come out the way "
        "they do — plus a downloadable PDF and Word copy."
    )
    tab_jump_button(
        "📄 Read the Full Paper — jump to that tab",
        key="jump_full_paper",
        target_tab_text="Full Paper",
    )

# ---------------------------------------------------------------------------
# Tab: full plain-language paper (read-in-app + downloadable PDF/Word)
# ---------------------------------------------------------------------------

with tab_paper:
    st.subheader("The Full Paper: A Plain-Language Guide")
    st.caption(
        "A complete, plain-language write-up of everything in this app — the "
        "biology, every method, every result, and why the numbers come out the way "
        "they do — written for a reader with no statistics or molecular-biology "
        "background, with every technical term defined before it's used. The "
        "executive summary below hits the highlights; download or open the full "
        "document (~7,000 words) for everything."
    )

    if full_paper_pdf_bytes is None and full_paper_docx_bytes is None:
        st.info(
            f"Paper files not found at `{FULL_PAPER_PDF_PATH}` / "
            f"`{FULL_PAPER_DOCX_PATH}`. Generate them from "
            f"`docs/full_methodology_paper.md` (see README's \"Running what exists\" "
            f"section)."
        )
    else:
        col1, col2, col3 = st.columns(3)
        if full_paper_pdf_bytes is not None:
            col1.download_button(
                "⬇️ Download as PDF",
                data=full_paper_pdf_bytes,
                file_name="neurospora_clock_network_full_paper.pdf",
                mime="application/pdf",
                width="stretch",
            )
        else:
            col1.info("PDF not generated yet.")

        if full_paper_docx_bytes is not None:
            col2.download_button(
                "⬇️ Download as Word (.docx)",
                data=full_paper_docx_bytes,
                file_name="neurospora_clock_network_full_paper.docx",
                mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                width="stretch",
            )
        else:
            col2.info("Word doc not generated yet.")

        if full_paper_pdf_bytes is not None:
            # A plain data: URI here gets silently blocked by Chromium's
            # top-level-navigation policy for data: URLs opened via
            # target="_blank" -- and a real server route would need Streamlit's
            # static-file serving enabled and a mirrored copy of the PDF kept in
            # sync, which is one more moving part than this needs. Instead:
            # embed the PDF bytes right in the page and, in the browser, decode
            # them into a Blob and open that via window.open() -- no server
            # route, no extra file, and blob: URLs aren't subject to the same
            # data: URI navigation block.
            #
            # This has to be rendered via components.html (a real iframe), not
            # st.markdown(unsafe_allow_html=True): Streamlit's markdown path
            # sanitizes out inline event-handler attributes like `onclick` even
            # with unsafe_allow_html=True (verified directly -- the attribute is
            # silently dropped), so the click handler would never fire there.
            # components.html's iframe isn't sanitized and its button click
            # keeps the real user-gesture context window.open() needs to avoid
            # being treated as a popup.
            #
            # Inside the onclick attribute specifically (not a <script> tag),
            # the browser wraps inline event-handler code in `with(document){
            # with(form){ with(element){ ... } } }` for legacy scoping reasons
            # -- and `document.URL` (a string, the page's own URL) shadows the
            # bare `URL` identifier that would otherwise resolve to the global
            # URL constructor. Referencing `window.URL.createObjectURL`
            # explicitly (verified directly -- the bare form throws
            # "URL.createObjectURL is not a function" here) sidesteps the
            # shadowing since `window` isn't a property of document/form/element.
            pdf_b64 = base64.b64encode(full_paper_pdf_bytes).decode("ascii")
            with col3:
                components.html(
                    f"""
                    <style>
                        /* The browser's UA stylesheet puts an 8px margin on body by
                        default -- left alone, it pushes the button down/right of
                        where the real st.download_button()s in the other two
                        columns sit, breaking the row's alignment. */
                        html, body {{ margin: 0; padding: 0; }}
                    </style>
                    <button onclick="
                        var bin = atob('{pdf_b64}');
                        var arr = new Uint8Array(bin.length);
                        for (var i = 0; i < bin.length; i++) {{ arr[i] = bin.charCodeAt(i); }}
                        var blob = new Blob([arr], {{type: 'application/pdf'}});
                        window.open(window.URL.createObjectURL(blob), '_blank');
                    " style="
                        display: inline-flex; align-items: center; justify-content: center;
                        width: 100%; box-sizing: border-box; padding: 0.4rem 0.75rem;
                        border: 1px solid rgba(49, 51, 63, 0.2); border-radius: 0.5rem;
                        background-color: #ffffff; color: inherit; text-decoration: none;
                        font-family: 'Source Sans', -apple-system, BlinkMacSystemFont,
                            'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
                        font-weight: 400; font-size: 16px; line-height: 1.6;
                        height: 2.5rem; cursor: pointer;
                    ">
                        🔗 Open PDF in New Tab
                    </button>
                    """,
                    height=40,
                )
        else:
            col3.info("PDF not generated yet.")

    st.divider()
    st.markdown("### Executive Summary")

    if full_paper_summary is None:
        st.info(f"Paper text not found at `{FULL_PAPER_MD_PATH}`.")
    else:
        st.markdown(full_paper_summary)
