#!/usr/bin/env python3

"""Build a combined bacterial + fungal antiSMASH overview page.

One row per BGC region, with summary statistics and charts on top and a
sortable/filterable table below. The page is self-contained (inline CSS, JS
and SVG) so it still works after being copied into the report's rsc/ tree.
"""

import html
import json
import re
from pathlib import Path
from urllib.parse import quote

import pandas as pd

# Same colours the Orochi report uses for its bacterial/fungal buttons.
COLOR_BACTERIA = "#1f77b4"
COLOR_FUNGI = "#03b48d"
COLOR_UNKNOWN = "#9e9e9e"

# Human-readable labels for the "taxon" column written by summarize_antismash.
TAXON_LABELS = {
    "bacteria": "Bacterial",
    "fungi": "Fungal",
}

# Longest bar chart; the remaining types are collapsed into one "other" bar.
MAX_TYPES_IN_CHART = 20


def safe_html(value):
    """Escape values before inserting them into HTML."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "N/A"

    text = str(value).strip()

    if not text or text.upper() == "NAN":
        return "N/A"

    return html.escape(text, quote=True)


def taxon_label(taxon):
    """Map the internal taxon key to the label shown in the page."""
    return TAXON_LABELS.get(str(taxon).strip().lower(), "Unknown")


def taxon_color(label):
    """Bar/segment colour for a taxon label."""
    if label == "Bacterial":
        return COLOR_BACTERIA

    if label == "Fungal":
        return COLOR_FUNGI

    return COLOR_UNKNOWN


def parse_record_data(text):
    """Extract the JSON array assigned to `var recordData` in regions.js."""
    match = re.search(r"var\s+recordData\s*=\s*", text)

    if not match:
        return None

    start = text.find("[", match.end())

    if start == -1:
        return None

    # regions.js holds several variables, so read to the matching bracket
    # instead of to the end of the file. Strings are tracked so that a "]"
    # inside a description does not end the array early.
    depth = 0
    in_string = False
    escaped = False

    for index in range(start, len(text)):
        character = text[index]

        if in_string:
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == '"':
                in_string = False

            continue

        if character == '"':
            in_string = True
        elif character == "[":
            depth += 1
        elif character == "]":
            depth -= 1

            if depth == 0:
                try:
                    return json.loads(text[start:index + 1])
                except json.JSONDecodeError:
                    return None

    return None


def load_anchors(regions_js_path):
    """Map (contig, region number) to the anchor used inside index.html.

    antiSMASH renders every region of a run into a single index.html and
    reaches them by fragment ("index.html#r1c1"). regions.js is what that
    page itself uses, so it is the reliable source for those anchors --
    record numbering cannot be recomputed safely from the summary table.
    """
    anchors = {}
    path = Path(regions_js_path)

    if not path.is_file():
        return anchors

    records = parse_record_data(
        path.read_text(encoding="utf-8", errors="replace")
    )

    if not isinstance(records, list):
        return anchors

    for record_index, record in enumerate(records, start=1):
        if not isinstance(record, dict):
            continue

        contig_keys = set()

        for field in ("seq_id", "id", "orig_id", "original_id"):
            value = record.get(field)

            if value:
                contig_keys.add(str(value).split()[0])

        regions = record.get("regions") or []

        for region_index, region in enumerate(regions, start=1):
            if not isinstance(region, dict):
                continue

            try:
                number = int(region.get("idx", region_index))
            except (TypeError, ValueError):
                number = region_index

            # Older antiSMASH versions omit "anchor"; its format is stable.
            anchor = region.get("anchor") or f"r{record_index}c{number}"

            for contig in contig_keys:
                anchors.setdefault((contig, number), str(anchor))

    return anchors


def load_table(table_path):
    """Read and normalise the combined BGC table."""
    try:
        table_df = pd.read_csv(table_path, sep="\t")
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Combined BGC table does not exist: {table_path}"
        ) from exc
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["bgc_id"])

    text_columns = [
        "bgc_id",
        "taxon",
        "contig_id",
        "bin_id",
        "taxonomy",
        "taxonomy_lowest",
        "bgc_product",
    ]

    for column in text_columns:
        if column not in table_df.columns:
            table_df[column] = "N/A"

        table_df[column] = (
            table_df[column]
            .fillna("N/A")
            .astype(str)
            .str.strip()
            .replace("", "N/A")
        )

    if table_df.empty:
        return table_df

    table_df["taxon_label"] = table_df["taxon"].apply(taxon_label)

    # Matches the key regions.js is indexed by (see load_anchors).
    table_df["contig_key"] = (
        table_df["contig_id"].str.split(n=1).str[0]
    )

    if "region_number" in table_df.columns:
        table_df["region_number"] = (
            pd.to_numeric(table_df["region_number"], errors="coerce")
            .fillna(1)
            .astype(int)
        )
    else:
        table_df["region_number"] = 1

    return table_df


def compute_stats(table_df):
    """Headline numbers shown in the cards above the charts."""
    if table_df.empty:
        return {
            "total": 0,
            "bacterial": 0,
            "fungal": 0,
            "binned": 0,
            "n_types": 0,
            "n_bins": 0,
        }

    binned = table_df[table_df["bin_id"] != "N/A"]

    return {
        "total": len(table_df),
        "bacterial": int((table_df["taxon_label"] == "Bacterial").sum()),
        "fungal": int((table_df["taxon_label"] == "Fungal").sum()),
        "binned": len(binned),
        "n_types": int(
            table_df.loc[table_df["bgc_product"] != "N/A", "bgc_product"]
            .nunique()
        ),
        "n_bins": int(binned["bin_id"].nunique()),
    }


def type_counts(table_df):
    """BGC counts per product type, split by taxon, most frequent first."""
    if table_df.empty:
        return []

    counts = (
        table_df.groupby(["bgc_product", "taxon_label"])
        .size()
        .unstack(fill_value=0)
    )

    counts["total"] = counts.sum(axis=1)
    counts = counts.sort_values("total", ascending=False)

    rows = []

    for product, row in counts.iterrows():
        rows.append(
            {
                "product": str(product),
                "total": int(row["total"]),
                "parts": [
                    (label, int(row[label]))
                    for label in ("Bacterial", "Fungal", "Unknown")
                    if label in counts.columns and int(row[label]) > 0
                ],
            }
        )

    return rows


def make_donut(stats):
    """Donut chart of bacterial vs fungal BGCs, drawn with dashed circles."""
    total = stats["total"]

    segments = [
        ("Bacterial", stats["bacterial"], COLOR_BACTERIA),
        ("Fungal", stats["fungal"], COLOR_FUNGI),
    ]

    other = total - stats["bacterial"] - stats["fungal"]

    if other > 0:
        segments.append(("Unknown", other, COLOR_UNKNOWN))

    radius = 60
    circumference = 2 * 3.141592653589793 * radius

    parts = []
    offset = 0.0

    for label, count, color in segments:
        if count <= 0:
            continue

        fraction = count / total if total else 0
        length = fraction * circumference

        parts.append(
            f'<circle class="donut-segment" cx="100" cy="100" r="{radius}" '
            f'fill="none" stroke="{color}" stroke-width="26" '
            f'stroke-dasharray="{length:.2f} {circumference - length:.2f}" '
            f'stroke-dashoffset="{-offset:.2f}">'
            f"<title>{safe_html(label)}: {count} BGCs</title>"
            "</circle>"
        )

        offset += length

    if not parts:
        parts.append(
            f'<circle cx="100" cy="100" r="{radius}" fill="none" '
            f'stroke="#e0e0e0" stroke-width="26"></circle>'
        )

    legend = []

    for label, count, color in segments:
        share = f"{100 * count / total:.0f}%" if total else "0%"

        legend.append(
            '<li>'
            f'<span class="swatch" style="background-color:{color}"></span>'
            f"{safe_html(label)}"
            f'<span class="legend-value">{count} ({share})</span>'
            "</li>"
        )

    return f"""
        <div class="chart-card">
            <h3>Bacterial vs fungal</h3>

            <div class="donut-wrap">
                <svg
                    viewBox="0 0 200 200"
                    class="donut"
                    role="img"
                    aria-label="Share of bacterial and fungal BGCs"
                >
                    <g transform="rotate(-90 100 100)">
                        {"".join(parts)}
                    </g>

                    <text
                        x="100"
                        y="94"
                        text-anchor="middle"
                        class="donut-total"
                    >{total}</text>

                    <text
                        x="100"
                        y="116"
                        text-anchor="middle"
                        class="donut-label"
                    >BGCs</text>
                </svg>

                <ul class="legend">
                    {"".join(legend)}
                </ul>
            </div>
        </div>
"""


def make_type_bars(rows):
    """Horizontal stacked bar chart of BGC counts per product type."""
    if not rows:
        return """
        <div class="chart-card">
            <h3>BGCs per type</h3>
            <p class="empty">No BGCs detected.</p>
        </div>
"""

    shown = rows[:MAX_TYPES_IN_CHART]
    hidden = rows[MAX_TYPES_IN_CHART:]

    if hidden:
        shown = shown + [
            {
                "product": f"other ({len(hidden)} types)",
                "total": sum(row["total"] for row in hidden),
                "parts": [
                    (
                        label,
                        sum(
                            count
                            for row in hidden
                            for part_label, count in row["parts"]
                            if part_label == label
                        ),
                    )
                    for label in ("Bacterial", "Fungal", "Unknown")
                ],
            }
        ]

        shown[-1]["parts"] = [
            (label, count)
            for label, count in shown[-1]["parts"]
            if count > 0
        ]

    # User units; the SVG itself scales to the card width via viewBox.
    label_width = 210
    bar_area = 520
    row_height = 24
    top_margin = 10

    height = top_margin * 2 + row_height * len(shown)
    max_total = max(row["total"] for row in shown) or 1

    elements = []

    for index, row in enumerate(shown):
        y = top_margin + index * row_height
        bar_y = y + 4
        bar_height = row_height - 10

        product = row["product"]
        display = product if len(product) <= 32 else product[:31] + "…"

        elements.append(
            f'<text x="{label_width - 8}" y="{y + row_height / 2 + 1}" '
            'text-anchor="end" class="bar-label">'
            f"{safe_html(display)}"
            f"<title>{safe_html(product)}</title>"
            "</text>"
        )

        x = label_width

        for label, count in row["parts"]:
            width = bar_area * count / max_total

            elements.append(
                f'<rect x="{x:.2f}" y="{bar_y}" width="{width:.2f}" '
                f'height="{bar_height}" fill="{taxon_color(label)}" '
                'rx="2">'
                f"<title>{safe_html(product)} &ndash; {safe_html(label)}: "
                f"{count}</title>"
                "</rect>"
            )

            x += width

        elements.append(
            f'<text x="{x + 6:.2f}" y="{y + row_height / 2 + 1}" '
            f'class="bar-value">{row["total"]}</text>'
        )

    legend = "".join(
        '<li>'
        f'<span class="swatch" style="background-color:{taxon_color(label)}">'
        "</span>"
        f"{label}"
        "</li>"
        for label in ("Bacterial", "Fungal")
    )

    return f"""
        <div class="chart-card chart-card-wide">
            <h3>BGCs per type</h3>

            <ul class="legend legend-inline">
                {legend}
            </ul>

            <svg
                viewBox="0 0 800 {height}"
                class="bars"
                role="img"
                aria-label="Number of BGCs per product type"
            >
                {"".join(elements)}
            </svg>
        </div>
"""


def make_select(column_index, label, values):
    """Dropdown filter for a low-cardinality column."""
    options = "".join(
        f'<option value="{safe_html(value)}">{safe_html(value)}</option>'
        for value in values
    )

    return f"""
                <th>
                    <select
                        class="column-filter"
                        data-column="{column_index}"
                        data-mode="exact"
                        aria-label="Filter by {safe_html(label)}"
                    >
                        <option value="">All {safe_html(label)}</option>
                        {options}
                    </select>
                </th>
"""


def make_text_filter(column_index, label, values, list_id):
    """Substring filter for a high-cardinality column, with suggestions."""
    options = "".join(
        f'<option value="{safe_html(value)}"></option>' for value in values
    )

    return f"""
                <th>
                    <input
                        type="search"
                        class="column-filter"
                        data-column="{column_index}"
                        data-mode="contains"
                        list="{list_id}"
                        placeholder="Filter {safe_html(label)}&hellip;"
                        aria-label="Filter by {safe_html(label)}"
                    >
                    <datalist id="{list_id}">{options}</datalist>
                </th>
"""


def unique_values(table_df, column):
    """Sorted unique non-missing values of a column, N/A kept last."""
    if table_df.empty or column not in table_df.columns:
        return []

    values = sorted(
        {value for value in table_df[column] if value and value != "N/A"},
        key=lambda value: value.lower(),
    )

    if (table_df[column] == "N/A").any():
        values.append("N/A")

    return values


def row_taxon(row):
    """Key of the antiSMASH run a BGC came from."""
    return "fungi" if row["taxon_label"] == "Fungal" else "bacteria"


def report_href(row, link_bases, anchors):
    """Link to the antiSMASH region page for one BGC, or None."""
    taxon = row_taxon(row)
    base = link_bases.get(taxon)

    if not base:
        return None

    anchor = anchors.get(taxon, {}).get(
        (row["contig_key"], int(row["region_number"]))
    )

    # Without a known anchor the report still opens, just at its first
    # region rather than at this one.
    fragment = f"#{quote(anchor)}" if anchor else ""

    return f"{base}/index.html{fragment}"


def per_bin_href(bin_id, link_bases):
    """Link to the regenerated per-bin antiSMASH report, or None."""
    base = link_bases.get("per_bin")

    if not base or bin_id == "N/A":
        return None

    return f"{base}/{quote(bin_id)}/index.html"


def link(href, label, title):
    """Anchor opening in a new tab, falling back to plain text."""
    if not href:
        return label

    return (
        f'<a href="{safe_html(href)}" target="_blank" '
        f'rel="noopener noreferrer" title="{safe_html(title)}">{label}</a>'
    )


def make_html(table_df, sample_pool, link_bases, anchors):
    """Generate the complete overview page."""
    sample_pool_html = safe_html(sample_pool)

    stats = compute_stats(table_df)
    rows = type_counts(table_df)

    cards = [
        ("Total BGCs", stats["total"]),
        ("Bacterial BGCs", stats["bacterial"]),
        ("Fungal BGCs", stats["fungal"]),
        ("BGCs in a MAG", stats["binned"]),
        ("MAGs with BGCs", stats["n_bins"]),
        ("Distinct BGC types", stats["n_types"]),
    ]

    cards_html = "".join(
        f"""
            <div class="stat-card">
                <div class="stat-value">{value}</div>
                <div class="stat-label">{safe_html(label)}</div>
            </div>
"""
        for label, value in cards
    )

    filters_html = (
        make_text_filter(0, "BGC id", unique_values(table_df, "bgc_id"), "filter-bgc-id")
        + make_select(1, "sources", unique_values(table_df, "taxon_label"))
        + make_text_filter(
            2, "taxonomy", unique_values(table_df, "taxonomy_lowest"), "filter-taxonomy"
        )
        + make_select(3, "bins", unique_values(table_df, "bin_id"))
        + make_select(4, "types", unique_values(table_df, "bgc_product"))
    )

    body_rows = []

    if table_df.empty:
        body_rows.append(
            """
            <tr>
                <td colspan="5">No BGCs detected for this sample pool</td>
            </tr>
"""
        )
    else:
        for _, row in table_df.iterrows():
            bgc_id = safe_html(row["bgc_id"])
            source = safe_html(row["taxon_label"])
            taxonomy_full = safe_html(row["taxonomy"])
            taxonomy_lowest = safe_html(row["taxonomy_lowest"])
            bin_id = safe_html(row["bin_id"])
            product = safe_html(row["bgc_product"])

            source_class = row["taxon_label"].lower()

            product_cell = (
                f'<span class="bgc-pill">{product}</span>'
                if product != "N/A"
                else "N/A"
            )

            bgc_cell = link(
                report_href(row, link_bases, anchors),
                bgc_id,
                f"Open this region in the {source.lower()} antiSMASH "
                f"report (contig {row['contig_id']})",
            )

            bin_cell = link(
                per_bin_href(row["bin_id"], link_bases),
                bin_id,
                f"Open the antiSMASH report for {row['bin_id']}",
            )

            body_rows.append(
                f"""
            <tr>
                <td data-value="{bgc_id}">{bgc_cell}</td>

                <td data-value="{source}">
                    <span class="tag tag-{source_class}">{source}</span>
                </td>

                <td data-value="{taxonomy_lowest}">
                    <em title="{taxonomy_full}">{taxonomy_lowest}</em>
                </td>

                <td data-value="{bin_id}">{bin_cell}</td>

                <td data-value="{product}">{product_cell}</td>
            </tr>
"""
            )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">

    <meta
        name="viewport"
        content="width=device-width, initial-scale=1.0"
    >

    <title>
        antiSMASH BGC overview - {sample_pool_html}
    </title>

    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            color: #222;
        }}

        h1 {{
            margin-bottom: 4px;
        }}

        h2 {{
            margin-top: 0;
            font-weight: normal;
            color: #555;
        }}

        h3 {{
            margin: 0 0 12px 0;
            font-size: 1em;
        }}

        .stats {{
            display: flex;
            flex-wrap: wrap;
            gap: 12px;
            margin: 20px 0;
        }}

        .stat-card {{
            flex: 1 1 140px;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 12px 16px;
            background-color: #fafafa;
        }}

        .stat-value {{
            font-size: 1.8em;
            font-weight: bold;
            color: #2f6f33;
        }}

        .stat-label {{
            font-size: 0.85em;
            color: #555;
        }}

        .charts {{
            display: flex;
            flex-wrap: wrap;
            gap: 12px;
            margin-bottom: 24px;
        }}

        .chart-card {{
            flex: 1 1 260px;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 16px;
        }}

        .chart-card-wide {{
            flex: 3 1 460px;
        }}

        .donut-wrap {{
            display: flex;
            align-items: center;
            gap: 16px;
            flex-wrap: wrap;
        }}

        .donut {{
            width: 170px;
            height: 170px;
            flex: 0 0 auto;
        }}

        .donut-total {{
            font-size: 26px;
            font-weight: bold;
            fill: #222;
        }}

        .donut-label {{
            font-size: 12px;
            fill: #666;
        }}

        .bars {{
            width: 100%;
            height: auto;
        }}

        .bar-label {{
            font-size: 12px;
            fill: #333;
        }}

        .bar-value {{
            font-size: 12px;
            fill: #666;
        }}

        .legend {{
            list-style: none;
            margin: 0;
            padding: 0;
            font-size: 0.9em;
        }}

        .legend li {{
            margin-bottom: 6px;
        }}

        .legend-inline {{
            display: flex;
            gap: 16px;
            margin-bottom: 8px;
        }}

        .legend-value {{
            color: #666;
            margin-left: 6px;
        }}

        .swatch {{
            display: inline-block;
            width: 12px;
            height: 12px;
            border-radius: 3px;
            margin-right: 6px;
            vertical-align: -1px;
        }}

        .toolbar {{
            display: flex;
            align-items: center;
            gap: 12px;
            flex-wrap: wrap;
            margin-bottom: 8px;
        }}

        #search {{
            flex: 1 1 280px;
            max-width: 420px;
            padding: 8px 10px;
            font-size: 1em;
            border: 1px solid #ccc;
            border-radius: 4px;
        }}

        #reset {{
            padding: 8px 14px;
            font-size: 0.9em;
            border: 1px solid #ccc;
            border-radius: 4px;
            background-color: #f5f5f5;
            cursor: pointer;
        }}

        #reset:hover {{
            background-color: #e9e9e9;
        }}

        #row-count {{
            color: #555;
            font-size: 0.9em;
        }}

        table {{
            border-collapse: collapse;
            width: 100%;
            margin-top: 8px;
        }}

        th,
        td {{
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
            vertical-align: top;
        }}

        thead tr.header-row th {{
            background-color: #4CAF50;
            color: white;
            cursor: pointer;
            user-select: none;
            position: relative;
            padding-right: 26px;
        }}

        thead tr.header-row th:hover {{
            background-color: #449d48;
        }}

        thead tr.header-row th.sort-asc::after {{
            content: "\\25B2";
            position: absolute;
            right: 8px;
        }}

        thead tr.header-row th.sort-desc::after {{
            content: "\\25BC";
            position: absolute;
            right: 8px;
        }}

        thead tr.filter-row th {{
            background-color: #f5f5f5;
            padding: 4px;
        }}

        .column-filter {{
            width: 100%;
            box-sizing: border-box;
            padding: 5px;
            font-size: 0.9em;
            border: 1px solid #ccc;
            border-radius: 3px;
            background-color: white;
        }}

        tbody tr:hover {{
            background-color: #f5f5f5;
        }}

        tbody a {{
            color: #0066cc;
            text-decoration: none;
        }}

        tbody a:hover {{
            text-decoration: underline;
        }}

        .bgc-pill {{
            display: inline-block;
            background-color: #e8f0fe;
            color: #1a3a6b;
            border-radius: 12px;
            padding: 2px 10px;
            font-size: 0.85em;
            white-space: nowrap;
        }}

        .tag {{
            display: inline-block;
            border-radius: 12px;
            padding: 2px 10px;
            font-size: 0.85em;
            color: white;
            white-space: nowrap;
        }}

        .tag-bacterial {{
            background-color: {COLOR_BACTERIA};
        }}

        .tag-fungal {{
            background-color: {COLOR_FUNGI};
        }}

        .tag-unknown {{
            background-color: {COLOR_UNKNOWN};
        }}

        .empty {{
            color: #777;
        }}
    </style>
</head>

<body>
    <h1>antiSMASH BGC overview</h1>
    <h2>Sample pool: {sample_pool_html}</h2>

    <div class="stats">
        {cards_html}
    </div>

    <div class="charts">
        {make_donut(stats)}
        {make_type_bars(rows)}
    </div>

    <div class="toolbar">
        <input
            type="search"
            id="search"
            placeholder="Search all columns&hellip;"
            aria-label="Search all columns"
        >

        <button type="button" id="reset">Reset filters</button>

        <span id="row-count"></span>
    </div>

    <table id="bgc-table">
        <thead>
            <tr class="header-row">
                <th data-type="text">BGC id</th>
                <th data-type="text">Source</th>
                <th data-type="text">Taxonomy</th>
                <th data-type="text">Bin</th>
                <th data-type="text">BGC type</th>
            </tr>

            <tr class="filter-row">
                {filters_html}
            </tr>
        </thead>

        <tbody>
            {"".join(body_rows)}
        </tbody>
    </table>

    <script>
        document.addEventListener("DOMContentLoaded", function () {{
            const table = document.getElementById("bgc-table");
            const tbody = table.querySelector("tbody");

            const headers = table.querySelectorAll(
                "thead tr.header-row th"
            );

            const filters = table.querySelectorAll(".column-filter");
            const search = document.getElementById("search");
            const resetButton = document.getElementById("reset");
            const rowCount = document.getElementById("row-count");

            const dataRows = Array.from(
                tbody.querySelectorAll("tr")
            ).filter(function (row) {{
                return !row.querySelector("td[colspan]");
            }});

            function cellValue(row, columnIndex) {{
                const cell = row.cells[columnIndex];

                if (!cell) {{
                    return "";
                }}

                const value = cell.dataset.value;

                return (
                    value !== undefined
                        ? value
                        : cell.textContent
                ).trim();
            }}

            function applyFilters() {{
                const searchTerm = search.value.trim().toLowerCase();

                let visible = 0;

                dataRows.forEach(function (row) {{
                    let keep = true;

                    filters.forEach(function (filter) {{
                        if (!keep) {{
                            return;
                        }}

                        const term = filter.value.trim().toLowerCase();

                        if (!term) {{
                            return;
                        }}

                        const value = cellValue(
                            row,
                            Number(filter.dataset.column)
                        ).toLowerCase();

                        if (filter.dataset.mode === "exact") {{
                            keep = value === term;
                        }} else {{
                            keep = value.indexOf(term) !== -1;
                        }}
                    }});

                    if (keep && searchTerm) {{
                        keep = row.textContent
                            .toLowerCase()
                            .indexOf(searchTerm) !== -1;
                    }}

                    row.style.display = keep ? "" : "none";

                    if (keep) {{
                        visible += 1;
                    }}
                }});

                rowCount.textContent =
                    "Showing " + visible + " of " +
                    dataRows.length + " BGCs";
            }}

            filters.forEach(function (filter) {{
                filter.addEventListener("input", applyFilters);
                filter.addEventListener("change", applyFilters);
            }});

            search.addEventListener("input", applyFilters);

            resetButton.addEventListener("click", function () {{
                filters.forEach(function (filter) {{
                    filter.value = "";
                }});

                search.value = "";

                applyFilters();
            }});

            headers.forEach(function (header, columnIndex) {{
                header.addEventListener("click", function () {{
                    if (dataRows.length === 0) {{
                        return;
                    }}

                    const dataType = header.dataset.type || "text";

                    const ascending =
                        !header.classList.contains("sort-asc");

                    headers.forEach(function (otherHeader) {{
                        otherHeader.classList.remove(
                            "sort-asc",
                            "sort-desc"
                        );
                    }});

                    header.classList.add(
                        ascending ? "sort-asc" : "sort-desc"
                    );

                    const sorted = dataRows.slice().sort(function (a, b) {{
                        const valueA = cellValue(a, columnIndex);
                        const valueB = cellValue(b, columnIndex);

                        const missingA =
                            valueA === "" ||
                            valueA.toUpperCase() === "N/A";

                        const missingB =
                            valueB === "" ||
                            valueB.toUpperCase() === "N/A";

                        // Keep missing values at the bottom.
                        if (missingA && missingB) {{
                            return 0;
                        }}

                        if (missingA) {{
                            return 1;
                        }}

                        if (missingB) {{
                            return -1;
                        }}

                        let comparison;

                        if (dataType === "number") {{
                            comparison = Number(valueA) - Number(valueB);
                        }} else {{
                            comparison = valueA.localeCompare(
                                valueB,
                                undefined,
                                {{
                                    numeric: true,
                                    sensitivity: "base"
                                }}
                            );
                        }}

                        return ascending ? comparison : -comparison;
                    }});

                    sorted.forEach(function (row) {{
                        tbody.appendChild(row);
                    }});
                }});
            }});

            applyFilters();
        }});
    </script>
</body>
</html>
"""


def main():
    table_path = Path(snakemake.input.combined_table)
    sample_pool = str(snakemake.params.sample_pool)

    antismash_dirs = dict(snakemake.params.antismash_dirs)
    link_bases = dict(snakemake.params.link_bases)

    table_df = load_table(table_path)

    # Kept per taxon: the two antiSMASH runs number their records
    # independently, so an anchor is only meaningful within its own report.
    anchors = {}

    for taxon, antismash_dir in antismash_dirs.items():
        anchors[taxon] = load_anchors(Path(antismash_dir) / "regions.js")

        if not anchors[taxon]:
            print(
                f"WARNING: no region anchors read from {antismash_dir}/"
                "regions.js -- BGC links for this taxon will open the "
                "report at its first region."
            )

    print(f"Sample pool: {sample_pool}")
    print(f"Number of BGC rows: {len(table_df)}")
    print(
        "Region anchors available: "
        + ", ".join(
            f"{taxon}={len(found)}" for taxon, found in anchors.items()
        )
    )

    if not table_df.empty:
        linked = sum(
            (row["contig_key"], int(row["region_number"]))
            in anchors.get(row_taxon(row), {})
            for _, row in table_df.iterrows()
        )

        print(f"BGCs linked to their own region: {linked}/{len(table_df)}")

    # The page is written twice: once next to the antiSMASH output and once
    # in the report's rsc/ tree, where the sibling report directories are
    # named differently. Only the link prefixes differ between the two.
    for output_name, bases in link_bases.items():
        output_path = Path(snakemake.output[output_name])

        html_content = make_html(
            table_df=table_df,
            sample_pool=sample_pool,
            link_bases=dict(bases),
            anchors=anchors,
        )

        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(html_content, encoding="utf-8")

        print(f"Generated combined BGC overview: {output_path}")


if __name__ == "__main__":
    main()