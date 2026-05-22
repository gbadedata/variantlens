"""
Plotly Dash Interactive Dashboard for VariantLens.

Panels:
  1. KPI cards       — total variants, pass rate, Ts/Tv, clinical count
  2. Consequence     — bar chart of variant consequence distribution
  3. AF distribution — allele frequency tier breakdown
  4. Variant class   — SNP vs indel composition
  5. Quality         — QUAL score distribution histogram
  6. Variant table   — filterable, sortable, paginated
  7. QC summary      — validation metrics panel
"""
import pandas as pd
import plotly.graph_objects as go
import dash
from dash import dcc, html, dash_table, Input, Output, callback
import dash_bootstrap_components as dbc

from src.processor import load_processed_data, load_summary

# ── Load data ────────────────────────────────────────────────────
df = load_processed_data()
summary = load_summary()

# ── Colours ──────────────────────────────────────────────────────
C = {
    "bg":        "#0D1117",
    "card":      "#161B22",
    "border":    "#30363D",
    "text":      "#E6EDF3",
    "subtext":   "#8B949E",
    "teal":      "#00D4AA",
    "teal_dark": "#008F75",
    "red":       "#F85149",
    "blue":      "#58A6FF",
    "yellow":    "#E3B341",
    "green":     "#3FB950",
    "purple":    "#BC8CFF",
    "orange":    "#FFA657",
}

CONSEQUENCE_COLOURS = {
    "synonymous":        C["green"],
    "missense":          C["orange"],
    "nonsense":          C["red"],
    "frameshift":        C["red"],
    "splice_region":     C["yellow"],
    "inframe_insertion": C["blue"],
    "inframe_deletion":  C["blue"],
    "MNP":               C["purple"],
    "complex_rearrangement": C["subtext"],
}

AF_COLOURS = {
    "rare":          C["blue"],
    "low_frequency": C["yellow"],
    "common":        C["red"],
    "unknown":       C["subtext"],
}

CLASS_COLOURS = {
    "SNP":       C["teal"],
    "insertion": C["blue"],
    "deletion":  C["red"],
    "MNP":       C["purple"],
    "complex":   C["subtext"],
}


# ── Helper components ────────────────────────────────────────────
def kpi_card(title, value, subtitle="", colour=None):
    colour = colour or C["teal"]
    return dbc.Card(
        dbc.CardBody([
            html.P(title, style={
                "color": C["subtext"], "fontSize": "0.72rem",
                "textTransform": "uppercase", "letterSpacing": "0.1em",
                "marginBottom": "4px",
            }),
            html.H3(value, style={
                "color": colour, "fontWeight": "700", "marginBottom": "2px",
            }),
            html.P(subtitle, style={
                "color": C["subtext"], "fontSize": "0.72rem", "marginBottom": 0,
            }),
        ]),
        style={
            "backgroundColor": C["card"],
            "border": f"1px solid {C['border']}",
            "borderRadius": "8px",
        },
    )


def qc_row(label, value):
    return html.Div([
        html.Span(label, style={"color": C["subtext"], "fontSize": "0.8rem"}),
        html.Span(value, style={
            "color": C["teal"], "fontSize": "0.8rem",
            "fontWeight": "600", "float": "right",
        }),
        html.Div(style={"clear": "both"}),
        html.Hr(style={"borderColor": C["border"], "margin": "5px 0"}),
    ])


def panel(title, children, height=None):
    style = {
        "backgroundColor": C["card"],
        "border": f"1px solid {C['border']}",
        "borderRadius": "8px",
        "marginBottom": "20px",
    }
    body_style = {}
    if height:
        body_style["height"] = height
    return dbc.Card([
        dbc.CardHeader(
            html.H6(title, style={"color": C["text"], "marginBottom": 0})
        ),
        dbc.CardBody(children, style=body_style),
    ], style=style)


# ── App ──────────────────────────────────────────────────────────
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.CYBORG],
    title="VariantLens",
)

# ── Layout ───────────────────────────────────────────────────────
app.layout = dbc.Container(fluid=True, style={
    "backgroundColor": C["bg"], "minHeight": "100vh", "padding": "24px",
}, children=[

    # Header
    dbc.Row([dbc.Col([
        html.H1("VariantLens", style={
            "color": C["teal"], "fontWeight": "700", "marginBottom": "4px",
        }),
        html.P(
            f"Genomic Variant Analysis · 1000 Genomes Project Phase 3 · "
            f"Chromosome 22 · {summary['reference_genome']} · "
            f"{summary['total_variants']} variants",
            style={"color": C["subtext"], "marginBottom": 0},
        ),
    ])], style={"marginBottom": "24px"}),

    # KPI cards
    dbc.Row([
        dbc.Col(kpi_card(
            "Total Variants", str(summary["total_variants"]),
            "passing QC validation"
        ), width=3),
        dbc.Col(kpi_card(
            "Pass Rate",
            f"{summary['validation']['pass_rate']}%",
            f"{summary['validation']['quarantined']} quarantined",
            colour=C["green"],
        ), width=3),
        dbc.Col(kpi_card(
            "Ts/Tv Ratio",
            str(summary["ts_tv_ratio"]),
            f"Ts={summary['ts_count']}  Tv={summary['tv_count']}",
            colour=C["blue"],
        ), width=3),
        dbc.Col(kpi_card(
            "SNPs",
            str(summary["variant_classes"].get("SNP", 0)),
            f"Indels: {summary.get('indel_count', 0)}",
            colour=C["teal"],
        ), width=3),
    ], style={"marginBottom": "24px"}),

    # Filters
    dbc.Row([dbc.Col([
        panel("Filters", dbc.Row([
            dbc.Col([
                html.Label("Variant Class:", style={
                    "color": C["subtext"], "fontSize": "0.85rem"
                }),
                dcc.Dropdown(
                    id="class-filter",
                    options=[{"label": "All classes", "value": "all"}] + [
                        {"label": c, "value": c}
                        for c in df["variant_class"].unique()
                    ],
                    value="all", clearable=False,
                    style={"backgroundColor": C["card"], "color": "#000"},
                ),
            ], width=3),
            dbc.Col([
                html.Label("Consequence:", style={
                    "color": C["subtext"], "fontSize": "0.85rem"
                }),
                dcc.Dropdown(
                    id="consequence-filter",
                    options=[{"label": "All consequences", "value": "all"}] + [
                        {"label": c, "value": c}
                        for c in sorted(df["consequence"].unique())
                    ],
                    value="all", clearable=False,
                    style={"backgroundColor": C["card"], "color": "#000"},
                ),
            ], width=3),
            dbc.Col([
                html.Label("AF Tier:", style={
                    "color": C["subtext"], "fontSize": "0.85rem"
                }),
                dcc.Dropdown(
                    id="af-tier-filter",
                    options=[{"label": "All tiers", "value": "all"}] + [
                        {"label": t, "value": t}
                        for t in ["rare", "low_frequency", "common"]
                    ],
                    value="all", clearable=False,
                    style={"backgroundColor": C["card"], "color": "#000"},
                ),
            ], width=3),
            dbc.Col([
                html.Label("PASS variants only:", style={
                    "color": C["subtext"], "fontSize": "0.85rem"
                }),
                dbc.Switch(id="pass-only-switch", value=True,
                           style={"marginTop": "8px"}),
            ], width=3),
        ])),
    ])], style={"marginBottom": "8px"}),

    # Charts row
    dbc.Row([
        # Consequence distribution
        dbc.Col([
            panel("Variant Consequences", [
                dcc.Graph(
                    id="consequence-chart",
                    style={"height": "280px"},
                    config={"displayModeBar": False},
                ),
            ]),
        ], width=5),

        # AF tier distribution
        dbc.Col([
            panel("Allele Frequency Tiers", [
                dcc.Graph(
                    id="af-tier-chart",
                    style={"height": "280px"},
                    config={"displayModeBar": False},
                ),
            ]),
        ], width=4),

        # QC summary
        dbc.Col([
            panel("QC Summary", [
                qc_row("Pass rate",
                       f"{summary['validation']['pass_rate']}%"),
                qc_row("Quarantined",
                       str(summary['validation']['quarantined'])),
                qc_row("Validation rules", "9 applied"),
                qc_row("Avg completeness",
                       f"{summary['validation']['avg_completeness']*100:.0f}%"),
                qc_row("Ts/Tv ratio", str(summary["ts_tv_ratio"])),
                qc_row("Reference", summary["reference_genome"]),
                qc_row("Source", "1000 Genomes Ph3"),
                qc_row("Chromosome", "22"),
            ]),
        ], width=3),
    ]),

    # Quality histogram
    dbc.Row([
        dbc.Col([
            panel("QUAL Score Distribution", [
                dcc.Graph(
                    id="qual-hist",
                    style={"height": "220px"},
                    config={"displayModeBar": False},
                ),
            ]),
        ], width=6),

        dbc.Col([
            panel("Variant Class Composition", [
                dcc.Graph(
                    id="class-donut",
                    style={"height": "220px"},
                    config={"displayModeBar": False},
                ),
            ]),
        ], width=6),
    ]),

    # Variant table
    dbc.Row([dbc.Col([
        panel("Variant Table", [
            html.Div(id="variant-table-container"),
        ]),
    ])]),

    # Footer
    dbc.Row([dbc.Col([
        html.Hr(style={"borderColor": C["border"]}),
        html.P(
            "VariantLens · github.com/gbadedata/variantlens · "
            "Python · Plotly Dash · FastAPI · pandas · "
            "Data: 1000 Genomes Project Phase 3",
            style={
                "color": C["subtext"], "fontSize": "0.75rem",
                "textAlign": "center",
            },
        ),
    ])], style={"marginTop": "16px"}),
])


# ── Callbacks ────────────────────────────────────────────────────

def apply_filters(class_f, cons_f, af_f, pass_only):
    filtered = df.copy()
    if class_f != "all":
        filtered = filtered[filtered["variant_class"] == class_f]
    if cons_f != "all":
        filtered = filtered[filtered["consequence"] == cons_f]
    if af_f != "all":
        filtered = filtered[filtered["af_tier"] == af_f]
    if pass_only:
        filtered = filtered[filtered["filter_outcome"] == "pass"]
    return filtered


@callback(
    Output("consequence-chart", "figure"),
    Input("class-filter", "value"),
    Input("af-tier-filter", "value"),
    Input("pass-only-switch", "value"),
)
def update_consequence_chart(class_f, af_f, pass_only):
    filtered = apply_filters(class_f, "all", af_f, pass_only)
    counts = filtered["consequence"].value_counts()
    colours = [
        CONSEQUENCE_COLOURS.get(c, C["subtext"]) for c in counts.index
    ]
    fig = go.Figure(go.Bar(
        x=list(counts.index),
        y=list(counts.values),
        marker_color=colours,
        hovertemplate="%{x}: %{y} variants<extra></extra>",
    ))
    fig.update_layout(
        plot_bgcolor=C["bg"], paper_bgcolor=C["card"],
        font=dict(color=C["text"], size=10),
        xaxis=dict(gridcolor=C["border"], tickangle=-30),
        yaxis=dict(gridcolor=C["border"], title="Count"),
        margin=dict(l=10, r=10, t=10, b=60),
        showlegend=False,
    )
    return fig


@callback(
    Output("af-tier-chart", "figure"),
    Input("class-filter", "value"),
    Input("consequence-filter", "value"),
    Input("pass-only-switch", "value"),
)
def update_af_chart(class_f, cons_f, pass_only):
    filtered = apply_filters(class_f, cons_f, "all", pass_only)
    tier_order = ["rare", "low_frequency", "common", "unknown"]
    counts = filtered["af_tier"].value_counts()
    tiers = [t for t in tier_order if t in counts.index]
    values = [int(counts[t]) for t in tiers]
    colours = [AF_COLOURS.get(t, C["subtext"]) for t in tiers]

    fig = go.Figure(go.Bar(
        x=tiers, y=values,
        marker_color=colours,
        hovertemplate="%{x}: %{y} variants<extra></extra>",
    ))
    fig.update_layout(
        plot_bgcolor=C["bg"], paper_bgcolor=C["card"],
        font=dict(color=C["text"], size=10),
        xaxis=dict(gridcolor=C["border"]),
        yaxis=dict(gridcolor=C["border"], title="Count"),
        margin=dict(l=10, r=10, t=10, b=40),
        showlegend=False,
    )
    return fig


@callback(
    Output("qual-hist", "figure"),
    Input("class-filter", "value"),
    Input("consequence-filter", "value"),
    Input("af-tier-filter", "value"),
    Input("pass-only-switch", "value"),
)
def update_qual_hist(class_f, cons_f, af_f, pass_only):
    filtered = apply_filters(class_f, cons_f, af_f, pass_only)
    qual_data = filtered["qual"].dropna()

    fig = go.Figure(go.Histogram(
        x=qual_data,
        nbinsx=40,
        marker_color=C["teal"],
        opacity=0.8,
        hovertemplate="QUAL %{x}: %{y} variants<extra></extra>",
    ))
    fig.update_layout(
        plot_bgcolor=C["bg"], paper_bgcolor=C["card"],
        font=dict(color=C["text"], size=10),
        xaxis=dict(gridcolor=C["border"], title="QUAL score"),
        yaxis=dict(gridcolor=C["border"], title="Count"),
        margin=dict(l=10, r=10, t=10, b=40),
        showlegend=False,
        bargap=0.05,
    )
    return fig


@callback(
    Output("class-donut", "figure"),
    Input("consequence-filter", "value"),
    Input("af-tier-filter", "value"),
    Input("pass-only-switch", "value"),
)
def update_class_donut(cons_f, af_f, pass_only):
    filtered = apply_filters("all", cons_f, af_f, pass_only)
    counts = filtered["variant_class"].value_counts()
    colours = [
        CLASS_COLOURS.get(c, C["subtext"]) for c in counts.index
    ]
    fig = go.Figure(go.Pie(
        labels=list(counts.index),
        values=list(counts.values),
        hole=0.55,
        marker=dict(colors=colours),
        hovertemplate="%{label}: %{value} (%{percent})<extra></extra>",
        textinfo="percent+label",
        textfont=dict(size=10, color=C["text"]),
    ))
    fig.update_layout(
        paper_bgcolor=C["card"],
        font=dict(color=C["text"]),
        margin=dict(l=10, r=10, t=10, b=10),
        showlegend=False,
    )
    return fig


@callback(
    Output("variant-table-container", "children"),
    Input("class-filter", "value"),
    Input("consequence-filter", "value"),
    Input("af-tier-filter", "value"),
    Input("pass-only-switch", "value"),
)
def update_table(class_f, cons_f, af_f, pass_only):
    filtered = apply_filters(class_f, cons_f, af_f, pass_only)
    cols = [
        "chrom", "pos", "variant_id", "ref", "alt",
        "qual", "af", "dp", "variant_class", "consequence",
        "gene", "af_tier", "ts_tv", "dbsnp_status",
    ]
    display = filtered[cols].copy()
    display["qual"] = display["qual"].round(1)
    display["af"] = display["af"].apply(
        lambda x: f"{x:.6f}" if x is not None and pd.notna(x) else ""
    )

    return dash_table.DataTable(
        data=display.to_dict("records"),
        columns=[
            {"name": c.replace("_", " ").title(), "id": c}
            for c in cols
        ],
        style_table={"overflowX": "auto"},
        style_cell={
            "backgroundColor": C["card"],
            "color": C["text"],
            "border": f"1px solid {C['border']}",
            "padding": "7px 10px",
            "fontSize": "0.82rem",
        },
        style_header={
            "backgroundColor": C["bg"],
            "color": C["teal"],
            "fontWeight": "600",
            "border": f"1px solid {C['border']}",
        },
        style_data_conditional=[
            {
                "if": {
                    "filter_query": '{consequence} = "missense"',
                    "column_id": "consequence",
                },
                "color": C["orange"],
                "fontWeight": "600",
            },
            {
                "if": {
                    "filter_query": '{consequence} = "nonsense"',
                    "column_id": "consequence",
                },
                "color": C["red"],
                "fontWeight": "600",
            },
            {
                "if": {
                    "filter_query": '{consequence} = "frameshift"',
                    "column_id": "consequence",
                },
                "color": C["red"],
                "fontWeight": "600",
            },
            {
                "if": {
                    "filter_query": '{af_tier} = "rare"',
                    "column_id": "af_tier",
                },
                "color": C["blue"],
            },
            {
                "if": {
                    "filter_query": '{af_tier} = "common"',
                    "column_id": "af_tier",
                },
                "color": C["red"],
            },
        ],
        page_size=15,
        sort_action="native",
        filter_action="native",
    )


if __name__ == "__main__":
    import structlog
    structlog.configure()
    app.run(debug=True, port=8050)
