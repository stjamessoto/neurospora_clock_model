"""Plotly visualizations for the clock network, for use in scripts and the
Streamlit app. networkx is used only for reference (layout is done by hand
below, mirroring the paper's Figure 1A hub-and-spoke diagram: WCC at the
center, the other 11 regulators arranged around it).

Colors follow a fixed, colorblind-validated palette (see dataviz skill /
references/palette.md): categorical hues assigned per individual regulator
(identity, REGULATOR_COLORS), the blue<->red diverging pair for
enriched/depleted (polarity), one sequential blue hue for magnitude/fit
lines, and neutral chart chrome (surface/gridline/ink) for everything that
isn't data.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from gene_layout_3d import GeneLayout3DResult
from network_builder import REGULATOR_ORDER, ClockNetwork
from network_stats import PAPER_POWER_LAW_EXPONENT, NetworkStatsResult
from null_models import NullModelResult

# --- palette tokens (light mode; see dataviz skill references/palette.md) ---
BLUE = "#2a78d6"        # categorical slot 1 / sequential hue / diverging pole "enriched"
RED = "#e34948"         # diverging pole "depleted"
SURFACE = "#fcfcfb"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"

# One distinct hue per regulator (11 total), extending the 8-slot categorical
# palette with 3 more hues chosen for CVD separation (validated with
# scripts/validate_palette.js --pairs all: worst all-pairs Machado-2009 CVD
# deltaE 11.2, the same floor the base 8 hues already sit at -- the 3
# additions don't make the worst case any worse). At 11 slots this sits in
# the palette's documented 8-12 "floor" band, legal only paired with
# secondary encoding -- every node this colors also carries a visible text
# label, in every figure that uses it.
REGULATOR_COLORS = dict(zip(REGULATOR_ORDER, [
    "#2a78d6",  # Wcc -- blue
    "#1baf7a",  # Reg1 NCU07392 (pro-1) -- aqua
    "#eda100",  # Reg2 NCU01640 (rpn-4) -- amber
    "#008300",  # Reg3 NCU06108 -- green
    "#4a3aa7",  # Reg4 NCU07155 -- violet
    "#e34948",  # PostReg5 NCU08295 (lhp-1) -- red
    "#e87ba4",  # PostReg6 NCU04504 (rrp-3) -- pink
    "#eb6834",  # PostReg7 NCU03363 -- orange
    "#72d025",  # PostReg8 NCU04799 (pab-1) -- lime
    "#39c6c6",  # PostReg9 NCU00919 (rok-1) -- cyan
    "#c825d0",  # PostReg10 NCU09349 (has-1) -- magenta-purple
]))
GENE_MULTI_COLOR = INK_MUTED  # genes with 2 regulators -- no single hub color applies

_AXIS_STYLE = dict(gridcolor=GRIDLINE, linecolor=BASELINE, tickfont=dict(color=INK_MUTED))


def hub_network_figure(
    net: ClockNetwork,
    enrichment_df: pd.DataFrame | None = None,
    q_threshold: float = 0.05,
    show_enriched: bool = True,
    show_depleted: bool = True,
    focus_label: str | None = None,
) -> go.Figure:
    """WCC-centered hub-and-spoke diagram (paper Fig. 1A layout).

    Neutral spokes = the fixed WCC -> regulator R structure (always
    present, reproduced from the paper). Diverging chords between the 10
    outer regulators = pairwise co-targeting from enrichment.py significant
    at ``q_threshold`` (a BH-adjusted q-value IS the smallest alpha at which
    a test rejects, so filtering on ``q_value <= q_threshold`` is exactly
    "significant at this FDR level" for any threshold, not just 0.05):
    blue = enriched (more shared targets than chance), red = depleted
    (fewer than chance) -- this project's novel extension. Node size scales
    with each regulator's target count (hub size); node color encodes
    regulator class (identity).

    ``focus_label`` (one of ``net.regulator_labels``, e.g. "PostReg5
    NCU08295") highlights that regulator and only its significant
    relationships at full opacity; everything else fades into the
    background rather than disappearing, so the surrounding structure is
    still visible for context.
    """
    labels = net.regulator_labels
    names = net.regulator_names
    classes = net.regulator_class
    hub_sizes = net.T.sum(axis=1)

    wcc_idx = labels.index("Wcc")
    others = [i for i in range(len(labels)) if i != wcc_idx]
    n_others = len(others)

    pos = {wcc_idx: (0.0, 0.0)}
    for k, i in enumerate(others):
        angle = 2 * np.pi * k / n_others - np.pi / 2
        pos[i] = (np.cos(angle), np.sin(angle))

    focus_is_set = focus_label is not None and focus_label in labels

    fig = go.Figure()

    # Fixed spokes: WCC -> every other regulator
    for i in others:
        x0, y0 = pos[wcc_idx]
        x1, y1 = pos[i]
        touches_focus = focus_is_set and labels[i] == focus_label
        opacity = 1.0 if (not focus_is_set or touches_focus) else 0.15
        fig.add_trace(go.Scatter(
            x=[x0, x1, None], y=[y0, y1, None], mode="lines",
            line=dict(color=BASELINE, width=1),
            opacity=opacity,
            hoverinfo="skip",
            showlegend=(i == others[0]),
            name="WCC regulation (fixed, from paper)",
            legendgroup="wcc-spoke",
        ))

    # Chords: pairwise enrichment/depletion significant at q_threshold (novel extension)
    if enrichment_df is not None and len(enrichment_df):
        label_to_idx = {lab: i for i, lab in enumerate(labels)}
        wanted_directions = set()
        if show_enriched:
            wanted_directions.add("enriched")
        if show_depleted:
            wanted_directions.add("depleted")
        pairs = enrichment_df[
            (enrichment_df["test_type"] == "pair")
            & (enrichment_df["q_value"] <= q_threshold)
            & (enrichment_df["direction"].isin(wanted_directions))
        ]
        legend_shown = {"enriched": False, "depleted": False}
        # draw dimmed non-focus chords first so focused ones land on top
        pairs = pairs.assign(
            _touches_focus=pairs.apply(
                lambda r: focus_is_set and (r["reg_a"] == focus_label or r["reg_b"] == focus_label), axis=1
            )
        ).sort_values("_touches_focus")
        for _, row in pairs.iterrows():
            a, b = row["reg_a"], row["reg_b"]
            if a not in label_to_idx or b not in label_to_idx:
                continue
            ia, ib = label_to_idx[a], label_to_idx[b]
            x0, y0 = pos[ia]
            x1, y1 = pos[ib]
            direction = row["direction"]
            color = BLUE if direction == "enriched" else RED
            q = max(row["q_value"], 1e-300)
            touches_focus = bool(row["_touches_focus"])
            if not focus_is_set or touches_focus:
                width = float(np.clip(-np.log10(q) / 8, 1, 6))
                opacity = 1.0
            else:
                width = 1.0
                opacity = 0.06
            show_legend = not legend_shown[direction]
            legend_shown[direction] = True
            fig.add_trace(go.Scatter(
                x=[x0, x1, None], y=[y0, y1, None], mode="lines",
                line=dict(color=color, width=width),
                opacity=opacity,
                hoverinfo="text",
                text=f"{row['regulators']}: {direction} (observed={row['observed_overlap']}, "
                     f"expected={row['expected_overlap']:.1f}, q={q:.2e})",
                name=f"significant: {direction}",
                showlegend=show_legend,
                legendgroup=direction,
            ))

    node_x = [pos[i][0] for i in range(len(labels))]
    node_y = [pos[i][1] for i in range(len(labels))]
    node_color = [REGULATOR_COLORS[lab] for lab in labels]
    max_hub = hub_sizes.max()
    node_size = [16 + 26 * (hub_sizes.iloc[i] / max_hub) for i in range(len(labels))]
    if focus_is_set:
        node_size = [
            s + 10 if labels[i] == focus_label else s
            for i, s in enumerate(node_size)
        ]
        node_opacity = [
            1.0 if (labels[i] == focus_label or labels[i] == "Wcc") else 0.3
            for i in range(len(labels))
        ]
        node_line_width = [3.0 if labels[i] == focus_label else 1.5 for i in range(len(labels))]
    else:
        node_opacity = [1.0] * len(labels)
        node_line_width = [1.5] * len(labels)
    hover = [
        f"{names[i]} ({labels[i]})<br>{int(hub_sizes.iloc[i])} targets<br>{classes[i]}"
        for i in range(len(labels))
    ]
    fig.add_trace(go.Scatter(
        x=node_x, y=node_y, mode="markers+text",
        marker=dict(
            size=node_size, color=node_color, opacity=node_opacity,
            line=dict(width=node_line_width, color=SURFACE),
        ),
        text=names, textposition="top center", textfont=dict(color=INK_PRIMARY),
        hoverinfo="text", hovertext=hover, showlegend=False,
    ))

    # dummy traces for a clean per-regulator legend
    for i, lab in enumerate(labels):
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers",
            marker=dict(size=10, color=REGULATOR_COLORS[lab]), name=names[i],
        ))

    title = "Clock network hub structure"
    if focus_is_set:
        idx = labels.index(focus_label)
        title += f" — focused on {names[idx]}"

    fig.update_layout(
        title=title,
        xaxis=dict(visible=False, range=[-1.4, 1.4]),
        yaxis=dict(visible=False, range=[-1.4, 1.4], scaleanchor="x"),
        height=650,
        plot_bgcolor=SURFACE,
        paper_bgcolor=SURFACE,
        font=dict(color=INK_SECONDARY),
        legend=dict(orientation="h", yanchor="bottom", y=-0.1),
    )
    return fig


def hub_network_figure_3d(
    net: ClockNetwork,
    enrichment_df: pd.DataFrame | None = None,
    q_threshold: float = 0.05,
    show_enriched: bool = True,
    show_depleted: bool = True,
    focus_label: str | None = None,
) -> go.Figure:
    """3D variant of hub_network_figure: WCC at the origin, the 5
    transcriptional regulators arranged in a ring above it and the 6
    post-transcriptional (RNA operon) regulators in a ring below it, so the
    two regulator classes separate visually along the vertical axis -- a
    genuinely 3D layout, not a 2D one just tilted. Drag to rotate, scroll to
    zoom. Same color/significance encoding and focus behavior as the 2D hub
    diagram; see that function's docstring for what the colors mean.

    Implementation note: unlike 2D Scatter, Plotly's Scatter3d only accepts
    a single scalar for marker.opacity and marker.line.width (no per-point
    arrays), so -- unlike the 2D version -- focus dimming is done by
    splitting nodes into a small number of traces (one per opacity/width
    combination) rather than one trace with per-point arrays.
    """
    labels = net.regulator_labels
    names = net.regulator_names
    classes = net.regulator_class
    hub_sizes = net.T.sum(axis=1)

    wcc_idx = labels.index("Wcc")
    others = [i for i in range(len(labels)) if i != wcc_idx]

    ring_groups: dict[str, list[int]] = {}
    for i in others:
        ring_groups.setdefault(classes[i], []).append(i)

    pos = {wcc_idx: (0.0, 0.0, 0.0)}
    for cls, z in [("transcriptional", 0.6), ("post-transcriptional", -0.6)]:
        idxs = ring_groups.get(cls, [])
        n = len(idxs)
        if n == 0:
            continue
        r = float(np.sqrt(max(1 - z ** 2, 0.05)))
        for k, i in enumerate(idxs):
            angle = 2 * np.pi * k / n
            pos[i] = (r * np.cos(angle), r * np.sin(angle), z)

    focus_is_set = focus_label is not None and focus_label in labels

    fig = go.Figure()

    # Fixed spokes: WCC -> every other regulator
    for i in others:
        x0, y0, z0 = pos[wcc_idx]
        x1, y1, z1 = pos[i]
        touches_focus = focus_is_set and labels[i] == focus_label
        opacity = 1.0 if (not focus_is_set or touches_focus) else 0.15
        fig.add_trace(go.Scatter3d(
            x=[x0, x1, None], y=[y0, y1, None], z=[z0, z1, None], mode="lines",
            line=dict(color=BASELINE, width=3),
            opacity=opacity,
            hoverinfo="skip",
            showlegend=(i == others[0]),
            name="WCC regulation (fixed, from paper)",
            legendgroup="wcc-spoke",
        ))

    # Chords: pairwise enrichment/depletion significant at q_threshold (novel extension)
    if enrichment_df is not None and len(enrichment_df):
        label_to_idx = {lab: i for i, lab in enumerate(labels)}
        wanted_directions = set()
        if show_enriched:
            wanted_directions.add("enriched")
        if show_depleted:
            wanted_directions.add("depleted")
        pairs = enrichment_df[
            (enrichment_df["test_type"] == "pair")
            & (enrichment_df["q_value"] <= q_threshold)
            & (enrichment_df["direction"].isin(wanted_directions))
        ]
        legend_shown = {"enriched": False, "depleted": False}
        pairs = pairs.assign(
            _touches_focus=pairs.apply(
                lambda r: focus_is_set and (r["reg_a"] == focus_label or r["reg_b"] == focus_label), axis=1
            )
        ).sort_values("_touches_focus")
        for _, row in pairs.iterrows():
            a, b = row["reg_a"], row["reg_b"]
            if a not in label_to_idx or b not in label_to_idx:
                continue
            ia, ib = label_to_idx[a], label_to_idx[b]
            x0, y0, z0 = pos[ia]
            x1, y1, z1 = pos[ib]
            direction = row["direction"]
            color = BLUE if direction == "enriched" else RED
            q = max(row["q_value"], 1e-300)
            touches_focus = bool(row["_touches_focus"])
            if not focus_is_set or touches_focus:
                width = float(np.clip(-np.log10(q) / 8, 1, 6)) + 1
                opacity = 1.0
            else:
                width = 2.0
                opacity = 0.06
            show_legend = not legend_shown[direction]
            legend_shown[direction] = True
            fig.add_trace(go.Scatter3d(
                x=[x0, x1, None], y=[y0, y1, None], z=[z0, z1, None], mode="lines",
                line=dict(color=color, width=width),
                opacity=opacity,
                hoverinfo="text",
                text=f"{row['regulators']}: {direction} (observed={row['observed_overlap']}, "
                     f"expected={row['expected_overlap']:.1f}, q={q:.2e})",
                name=f"significant: {direction}",
                showlegend=show_legend,
                legendgroup=direction,
            ))

    # Nodes, grouped into a handful of traces so opacity/line-width (scalar-only
    # in Scatter3d) can still vary between the focused node, WCC, and the rest.
    max_hub = hub_sizes.max()

    def _node_size(i: int) -> float:
        base = 10 + 20 * (hub_sizes.iloc[i] / max_hub)
        return base + 8 if focus_is_set and labels[i] == focus_label else base

    if focus_is_set:
        focus_idxs = [i for i in range(len(labels)) if labels[i] == focus_label]
        wcc_idxs = [i for i in range(len(labels)) if labels[i] == "Wcc" and labels[i] != focus_label]
        dim_idxs = [i for i in range(len(labels)) if i not in focus_idxs and i not in wcc_idxs]
        node_groups = [
            (focus_idxs, 1.0, 3.0),
            (wcc_idxs, 1.0, 1.5),
            (dim_idxs, 0.25, 1.0),
        ]
    else:
        node_groups = [([i for i in range(len(labels))], 1.0, 1.5)]

    for idxs, opacity, line_width in node_groups:
        if not idxs:
            continue
        fig.add_trace(go.Scatter3d(
            x=[pos[i][0] for i in idxs], y=[pos[i][1] for i in idxs], z=[pos[i][2] for i in idxs],
            mode="markers+text",
            marker=dict(
                size=[_node_size(i) for i in idxs],
                color=[REGULATOR_COLORS[labels[i]] for i in idxs],
                opacity=opacity,
                line=dict(width=line_width, color=SURFACE),
            ),
            text=[names[i] for i in idxs], textposition="top center",
            textfont=dict(color=INK_PRIMARY, size=11),
            hoverinfo="text",
            hovertext=[
                f"{names[i]} ({labels[i]})<br>{int(hub_sizes.iloc[i])} targets<br>{classes[i]}"
                for i in idxs
            ],
            showlegend=False,
        ))

    # dummy traces for a clean per-regulator legend
    for i, lab in enumerate(labels):
        fig.add_trace(go.Scatter3d(
            x=[None], y=[None], z=[None], mode="markers",
            marker=dict(size=6, color=REGULATOR_COLORS[lab]), name=names[i],
        ))

    title = "Clock network hub structure (3D — drag to rotate, scroll to zoom)"
    if focus_is_set:
        idx = labels.index(focus_label)
        title += f" — focused on {names[idx]}"

    axis_kwargs = dict(visible=False, showbackground=False, showgrid=False, zeroline=False, range=[-1.3, 1.3])
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis=axis_kwargs, yaxis=axis_kwargs, zaxis=axis_kwargs,
            camera=dict(eye=dict(x=1.5, y=1.5, z=0.9)),
            bgcolor=SURFACE,
        ),
        height=650,
        paper_bgcolor=SURFACE,
        font=dict(color=INK_SECONDARY),
        legend=dict(orientation="h", yanchor="bottom", y=-0.05),
        margin=dict(l=0, r=0, t=40, b=0),
    )
    return fig


def ffl_network_figure(net: ClockNetwork, ffl_result, focus_label: str | None = None) -> go.Figure:
    """2D node-link diagram of every feedforward loop counted by
    ``ffl_analysis.count_ffls``: WCC (A) at the center, every regulator (B)
    that completes at least one FFL arranged in a ring around it (gray
    spokes, the fixed WCC -> B structure), and each B's FFL-completing
    target genes (C) as small dots hanging off it -- solid line for the
    B -> C edge, dashed line for the WCC -> C edge that's the second leg of
    the same triangle. This doesn't introduce any new data: it's the same
    ``(R @ T) * T`` result already shown as numbers in the FFLs & Null Model
    tab, laid out as the loops it's actually counting.
    """
    labels = net.regulator_labels
    names = net.regulator_names
    classes = net.regulator_class
    wcc_idx = labels.index("Wcc")

    # ffl_result.ffl_matrix is (R @ T) * T -- since R's only nonzero row is
    # WCC's, every nonzero cell of that matrix sits on the WCC row (see
    # ffl_analysis.py's docstring), not on each B regulator's own row. The
    # actual per-(B, C) FFL pairs -- gene c completes a loop through
    # regulator b iff BOTH WCC and b target c -- are recovered the same way
    # ffl_analysis.py computes per_regulator_ffl_count: an AND of the two
    # regulators' rows in net.T, not a lookup into ffl_matrix by row.
    wcc_row = net.T.loc["Wcc"].astype(bool)
    b_labels = [lab for lab in labels if lab != "Wcc" and (wcc_row & net.T.loc[lab].astype(bool)).any()]
    n_b = len(b_labels)

    b_pos = {}
    for k, lab in enumerate(b_labels):
        angle = 2 * np.pi * k / max(n_b, 1) - np.pi / 2
        b_pos[lab] = (2.2 * np.cos(angle), 2.2 * np.sin(angle))
    wcc_pos = (0.0, 0.0)

    gene_pos = {}
    gene_of = {}  # gene -> b_label
    for lab in b_labels:
        genes = list(net.T.columns[wcc_row & net.T.loc[lab].astype(bool)])
        n = len(genes)
        bx, by = b_pos[lab]
        angle0 = np.arctan2(by, bx)
        spread = np.linspace(-0.35, 0.35, n) if n > 1 else [0.0]
        for gene, da in zip(genes, spread):
            a = angle0 + da
            gene_pos[gene] = (3.3 * np.cos(a), 3.3 * np.sin(a))
            gene_of[gene] = lab

    focus_is_set = focus_label is not None and focus_label in labels
    fig = go.Figure()

    # WCC -> B spokes (fixed structure)
    for lab in b_labels:
        x1, y1 = b_pos[lab]
        touches_focus = focus_is_set and lab == focus_label
        opacity = 1.0 if (not focus_is_set or touches_focus) else 0.15
        fig.add_trace(go.Scatter(
            x=[wcc_pos[0], x1, None], y=[wcc_pos[1], y1, None], mode="lines",
            line=dict(color=BASELINE, width=2),
            opacity=opacity, hoverinfo="skip", showlegend=(lab == b_labels[0]),
            name="WCC → B (fixed)", legendgroup="wcc-b",
        ))

    # B -> C (solid) and WCC -> C (dashed) -- the two legs completing each loop
    legend_shown = {"b_to_c": False, "wcc_to_c": False}
    for gene, lab in gene_of.items():
        gx, gy = gene_pos[gene]
        bx, by = b_pos[lab]
        touches_focus = focus_is_set and lab == focus_label
        opacity = 1.0 if (not focus_is_set or touches_focus) else 0.06
        color = REGULATOR_COLORS[lab]
        fig.add_trace(go.Scatter(
            x=[bx, gx, None], y=[by, gy, None], mode="lines",
            line=dict(color=color, width=1.5),
            opacity=opacity, hoverinfo="text", text=f"{lab} → {gene}",
            showlegend=not legend_shown["b_to_c"], name="B → C", legendgroup="b_to_c",
        ))
        legend_shown["b_to_c"] = True
        fig.add_trace(go.Scatter(
            x=[wcc_pos[0], gx, None], y=[wcc_pos[1], gy, None], mode="lines",
            line=dict(color=REGULATOR_COLORS["Wcc"], width=1, dash="dot"),
            opacity=opacity, hoverinfo="text", text=f"Wcc → {gene}",
            showlegend=not legend_shown["wcc_to_c"], name="WCC → C", legendgroup="wcc_to_c",
        ))
        legend_shown["wcc_to_c"] = True

    # gene markers
    for lab in b_labels:
        genes = [g for g, b in gene_of.items() if b == lab]
        if not genes:
            continue
        touches_focus = focus_is_set and lab == focus_label
        opacity = 1.0 if (not focus_is_set or touches_focus) else 0.1
        color = REGULATOR_COLORS[lab]
        fig.add_trace(go.Scatter(
            x=[gene_pos[g][0] for g in genes], y=[gene_pos[g][1] for g in genes],
            mode="markers", marker=dict(size=8, color=color, opacity=opacity, line=dict(width=1, color=SURFACE)),
            hoverinfo="text", hovertext=genes, showlegend=False,
        ))

    # WCC + B regulator nodes
    node_labels = ["Wcc"] + b_labels
    node_pos = {**{"Wcc": wcc_pos}, **b_pos}
    hub_sizes = net.T.sum(axis=1)
    max_hub = hub_sizes.max()
    node_x = [node_pos[lab][0] for lab in node_labels]
    node_y = [node_pos[lab][1] for lab in node_labels]
    node_color = [REGULATOR_COLORS[lab] for lab in node_labels]
    node_size = [16 + 26 * (hub_sizes.loc[lab] / max_hub) for lab in node_labels]
    if focus_is_set:
        node_size = [s + 10 if lab == focus_label else s for s, lab in zip(node_size, node_labels)]
        node_opacity = [1.0 if (lab == focus_label or lab == "Wcc") else 0.3 for lab in node_labels]
    else:
        node_opacity = [1.0] * len(node_labels)
    fig.add_trace(go.Scatter(
        x=node_x, y=node_y, mode="markers+text",
        marker=dict(size=node_size, color=node_color, opacity=node_opacity, line=dict(width=1.5, color=SURFACE)),
        text=[names[labels.index(lab)] for lab in node_labels], textposition="top center",
        textfont=dict(color=INK_PRIMARY),
        hoverinfo="text",
        hovertext=[
            f"{names[labels.index(lab)]} ({lab})<br>"
            f"{int(ffl_result.per_regulator_ffl_count.get(lab, 0))} FFLs as B"
            if lab != "Wcc" else f"WCC — {ffl_result.total_ffl_count} total FFLs (always A)"
            for lab in node_labels
        ],
        showlegend=False,
    ))

    title = f"Feedforward-loop network ({ffl_result.total_ffl_count} FFLs, {len(gene_pos)} target genes)"
    if focus_is_set:
        title += f" — focused on {names[labels.index(focus_label)]}"

    fig.update_layout(
        title=title,
        xaxis=dict(visible=False, range=[-4, 4]),
        yaxis=dict(visible=False, range=[-4, 4], scaleanchor="x"),
        height=650,
        plot_bgcolor=SURFACE,
        paper_bgcolor=SURFACE,
        font=dict(color=INK_SECONDARY),
        legend=dict(orientation="h", yanchor="bottom", y=-0.1),
    )
    return fig


def regulator_gene_figure_3d(
    net: ClockNetwork,
    layout: GeneLayout3DResult,
    focus_label: str | None = None,
) -> go.Figure:
    """3D counterpart to the 2D Cytoscape.js force-directed Regulators +
    Genes graph in the same subtab: same graph, same per-regulator colors
    (``REGULATOR_COLORS``), but positions come from
    ``gene_layout_3d.compute_layout`` -- a real 3D force-directed layout
    (networkx's ``spring_layout(dim=3)``, precomputed and cached, ~45s to
    (re)compute) rather than Cytoscape's client-side "fcose". The physics
    itself pulls each regulator's exclusive genes into their own cluster
    purely from shared-edge attraction, the same way the 2D view does.

    ``focus_label`` isolates one regulator's genes at full opacity and
    fades the rest -- Plotly has no built-in click-to-focus the way the
    Cytoscape component does, so this is a manual dropdown instead.
    """
    pos_by_id = {nid: layout.positions[i] for i, nid in enumerate(layout.node_ids)}
    labels = net.regulator_labels
    names = net.regulator_names
    hub_sizes = net.T.sum(axis=1)

    reg_pos = {lab: pos_by_id[f"reg::{lab}"] for lab in labels}
    gene_regs = {gene: list(net.T.index[net.T[gene].astype(bool)]) for gene in net.T.columns}

    focus_is_set = focus_label is not None and focus_label in labels
    fig = go.Figure()

    # regulator -> gene edges, one trace per regulator for focus dimming
    for lab in labels:
        genes = [g for g, regs in gene_regs.items() if lab in regs]
        if not genes:
            continue
        rx, ry, rz = reg_pos[lab]
        xs, ys, zs = [], [], []
        for gene in genes:
            gx, gy, gz = pos_by_id[f"gene::{gene}"]
            xs += [rx, gx, None]
            ys += [ry, gy, None]
            zs += [rz, gz, None]
        touches_focus = focus_is_set and lab == focus_label
        opacity = (1.0 if touches_focus else 0.03) if focus_is_set else 0.15
        fig.add_trace(go.Scatter3d(
            x=xs, y=ys, z=zs, mode="lines",
            line=dict(color=REGULATOR_COLORS[lab], width=1),
            opacity=opacity, hoverinfo="skip", showlegend=False,
        ))

    # gene markers, grouped by single regulator (colored to match) vs. multi (neutral)
    marker_groups: dict[str, list[str]] = {}
    for gene, regs in gene_regs.items():
        key = regs[0] if len(regs) == 1 else "__multi__"
        marker_groups.setdefault(key, []).append(gene)

    for key, genes in marker_groups.items():
        color = GENE_MULTI_COLOR if key == "__multi__" else REGULATOR_COLORS[key]
        touches_focus = focus_is_set and key == focus_label
        opacity = (0.85 if touches_focus else 0.05) if focus_is_set else 0.55
        xs = [pos_by_id[f"gene::{g}"][0] for g in genes]
        ys = [pos_by_id[f"gene::{g}"][1] for g in genes]
        zs = [pos_by_id[f"gene::{g}"][2] for g in genes]
        fig.add_trace(go.Scatter3d(
            x=xs, y=ys, z=zs, mode="markers",
            marker=dict(size=2.5, color=color, opacity=opacity),
            hoverinfo="text",
            hovertext=[f"{g}<br>regulated by: {', '.join(gene_regs[g])}" for g in genes],
            showlegend=False,
        ))

    # regulator nodes
    max_hub = hub_sizes.max()

    def node_size(lab: str) -> float:
        base = 10 + 20 * (hub_sizes.loc[lab] / max_hub)
        return base + 8 if focus_is_set and lab == focus_label else base

    if focus_is_set:
        focus_labs = [focus_label]
        wcc_labs = [l for l in labels if l == "Wcc" and l != focus_label]
        dim_labs = [l for l in labels if l not in focus_labs and l not in wcc_labs]
        node_groups = [(focus_labs, 1.0, 3.0), (wcc_labs, 1.0, 1.5), (dim_labs, 0.3, 1.0)]
    else:
        node_groups = [(labels, 1.0, 1.5)]

    for labs, opacity, line_width in node_groups:
        if not labs:
            continue
        fig.add_trace(go.Scatter3d(
            x=[reg_pos[l][0] for l in labs], y=[reg_pos[l][1] for l in labs], z=[reg_pos[l][2] for l in labs],
            mode="markers+text",
            marker=dict(
                size=[node_size(l) for l in labs],
                color=[REGULATOR_COLORS[l] for l in labs],
                opacity=opacity,
                line=dict(width=line_width, color=SURFACE),
            ),
            text=[names[labels.index(l)] for l in labs], textposition="top center",
            textfont=dict(color=INK_PRIMARY, size=11),
            hoverinfo="text",
            hovertext=[
                f"{names[labels.index(l)]} ({l})<br>{int(hub_sizes.loc[l])} targets" for l in labs
            ],
            showlegend=False,
        ))

    # dummy traces for a clean per-regulator legend
    for lab in labels:
        fig.add_trace(go.Scatter3d(
            x=[None], y=[None], z=[None], mode="markers",
            marker=dict(size=6, color=REGULATOR_COLORS[lab]), name=names[labels.index(lab)],
        ))

    title = f"Regulators + all {len(gene_regs):,} target genes (3D force-directed — drag to rotate, scroll to zoom)"
    if focus_is_set:
        title += f" — focused on {names[labels.index(focus_label)]}"

    axis_kwargs = dict(visible=False, showbackground=False, showgrid=False, zeroline=False)
    fig.update_layout(
        title=title,
        scene=dict(xaxis=axis_kwargs, yaxis=axis_kwargs, zaxis=axis_kwargs, bgcolor=SURFACE),
        height=700,
        paper_bgcolor=SURFACE,
        font=dict(color=INK_SECONDARY),
        legend=dict(orientation="h", yanchor="bottom", y=-0.05),
        margin=dict(l=0, r=0, t=40, b=0),
    )
    return fig


def degree_distribution_figure(result: NetworkStatsResult) -> go.Figure:
    """Log-log degree-value frequency plot with the fitted power law overlaid.

    Compare the fitted slope to the paper's reported -1.4 (Section V-E). The
    observed points and their own fitted trend share the sequential blue
    hue (same entity); the paper's reference exponent is neutral gray,
    since it isn't this project's data.
    """
    dd = result.degree_distribution
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=dd.index, y=dd.values, mode="markers", name="observed (this export)",
        marker=dict(size=9, color=BLUE),
    ))

    x_fit = np.array([dd.index.min(), dd.index.max()], dtype=float)
    y_fit = (10 ** result.power_law_intercept) * x_fit ** result.power_law_exponent
    fig.add_trace(go.Scatter(
        x=x_fit, y=y_fit, mode="lines",
        line=dict(color=BLUE, dash="solid", width=2),
        name=f"fitted exponent = {result.power_law_exponent:.2f} (R2={result.power_law_r_squared:.2f})",
    ))

    y_paper = (10 ** result.power_law_intercept) * x_fit ** PAPER_POWER_LAW_EXPONENT
    fig.add_trace(go.Scatter(
        x=x_fit, y=y_paper, mode="lines",
        line=dict(color=INK_MUTED, dash="dash", width=2),
        name=f"paper's reported exponent = {PAPER_POWER_LAW_EXPONENT}",
    ))

    fig.update_layout(
        title="Node degree distribution (regulators + target genes), log-log",
        xaxis=dict(title="degree (connections per node)", type="log", **_AXIS_STYLE),
        yaxis=dict(title="number of nodes", type="log", **_AXIS_STYLE),
        height=500,
        plot_bgcolor=SURFACE,
        paper_bgcolor=SURFACE,
        font=dict(color=INK_SECONDARY),
    )
    return fig


def null_distribution_figure(null_result: NullModelResult) -> go.Figure:
    """Histogram of FFL counts under degree-preserving rewiring, observed value marked.

    The null distribution itself uses the sequential blue hue (magnitude);
    the observed value is a single pinned reference, marked in red to draw
    the eye, matching the "depleted" pole used elsewhere for a below-chance
    reading (the observed count sits inside, not beyond, this null body).
    """
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=null_result.null_counts, nbinsx=30, name="null distribution (rewired networks)",
        marker_color=BLUE, opacity=0.75,
    ))
    fig.add_vline(
        x=null_result.observed_ffl_count, line_dash="dash", line_color=RED, line_width=2,
        annotation_text=f"observed = {null_result.observed_ffl_count}",
        annotation_position="top right", annotation_font_color=INK_PRIMARY,
    )
    fig.update_layout(
        title=(
            f"FFL count: observed vs. degree-preserving null "
            f"(Z={null_result.z_score:.2f}, two-sided p={null_result.p_value_two_sided:.3g})"
        ),
        xaxis=dict(title="FFL count in randomized network", **_AXIS_STYLE),
        yaxis=dict(title="frequency", **_AXIS_STYLE),
        height=450,
        plot_bgcolor=SURFACE,
        paper_bgcolor=SURFACE,
        font=dict(color=INK_SECONDARY),
    )
    return fig


def binding_strength_heatmap(net: ClockNetwork, table: pd.DataFrame) -> go.Figure:
    """Heatmap of the function x regulator binding-strength table (joint
    VTENS/MCMC point estimates, from ``binding_strength.py`` -- see that
    module's docstring for why these are joint estimates, not independent
    per-parameter distributions, and shouldn't be read as the result of a
    "does the ensemble exclude zero" test).

    Sequential blue (magnitude) encodes strength; a cell of exactly 0.0
    ("not the estimated dominant regulator for this function") renders as
    blank chart surface rather than a pale blue, so true zeros read as an
    absence, not a low value -- there's a real gap in this table between 0
    and the smallest observed positive strength (~0.36), so this is a
    faithful cutoff, not a rounding trick.
    """
    z = table.values
    zmax = float(z.max())
    eps = 0.06  # position fraction below which cells render as blank surface

    colorscale = [
        [0.0, SURFACE],
        [eps, SURFACE],
        [eps, "#cde2fb"],
        [0.4, "#5598e7"],
        [0.7, BLUE],
        [1.0, "#0d366b"],
    ]

    fig = go.Figure(go.Heatmap(
        z=z,
        x=list(net.regulator_names),
        y=list(table.index),
        colorscale=colorscale,
        zmin=0, zmax=zmax,
        xgap=2, ygap=2,
        colorbar=dict(title="binding<br>strength", tickfont=dict(color=INK_MUTED)),
        hovertemplate="%{y} × %{x}<br>binding strength: %{z:.3f}<extra></extra>",
    ))
    fig.update_layout(
        title="Binding strength by function (blank = not the estimated dominant regulator)",
        xaxis=dict(side="top", tickangle=-45, **_AXIS_STYLE),
        yaxis=dict(autorange="reversed", **_AXIS_STYLE),
        height=650,
        plot_bgcolor=SURFACE,
        paper_bgcolor=SURFACE,
        font=dict(color=INK_SECONDARY),
        margin=dict(l=220, t=140),
    )
    return fig
