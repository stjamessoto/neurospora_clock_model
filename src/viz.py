"""Plotly visualizations for the clock network, for use in scripts and the
Streamlit app. networkx is used only for reference (layout is done by hand
below, mirroring the paper's Figure 1A hub-and-spoke diagram: WCC at the
center, the other 11 regulators arranged around it).

Colors follow a fixed, colorblind-validated palette (see dataviz skill /
references/palette.md): categorical hues assigned by regulator class
(identity), the blue<->red diverging pair for enriched/depleted (polarity),
one sequential blue hue for magnitude/fit lines, and neutral chart chrome
(surface/gridline/ink) for everything that isn't data.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from network_builder import ClockNetwork
from network_stats import PAPER_POWER_LAW_EXPONENT, NetworkStatsResult
from null_models import NullModelResult

# --- palette tokens (light mode; see dataviz skill references/palette.md) ---
BLUE = "#2a78d6"        # categorical slot 1 / sequential hue / diverging pole "enriched"
GREEN = "#008300"       # categorical slot 4 -- matches paper's Fig. 1A transcriptional-regulator green
VIOLET = "#4a3aa7"      # categorical slot 5 -- master regulator (WCC), distinct from the two classes
ORANGE = "#eb6834"      # categorical slot 8 -- matches paper's Fig. 1A post-transcriptional-regulator orange
RED = "#e34948"         # diverging pole "depleted"
SURFACE = "#fcfcfb"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"

CLASS_COLORS = {
    "master": VIOLET,
    "transcriptional": GREEN,
    "post-transcriptional": ORANGE,
}
CLASS_LABELS = {
    "master": "Master regulator (WCC)",
    "transcriptional": "Transcriptional regulator",
    "post-transcriptional": "Post-transcriptional regulator (RNA operon)",
}

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
    node_color = [CLASS_COLORS[c] for c in classes]
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

    # dummy traces for a clean class legend
    for cls, color in CLASS_COLORS.items():
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers",
            marker=dict(size=10, color=color), name=CLASS_LABELS[cls],
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
                color=[CLASS_COLORS[classes[i]] for i in idxs],
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

    # dummy traces for a clean class legend
    for cls, color in CLASS_COLORS.items():
        fig.add_trace(go.Scatter3d(
            x=[None], y=[None], z=[None], mode="markers",
            marker=dict(size=6, color=color), name=CLASS_LABELS[cls],
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
