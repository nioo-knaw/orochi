#!/usr/bin/env python3

import html
from pathlib import Path

import pandas as pd


def safe_html(value):
    """Escape values before inserting them into HTML."""
    if pd.isna(value):
        return "N/A"

    return html.escape(str(value), quote=True)


def bgc_types_to_pills(bgc_types_str):
    """
    Convert a string such as:

        NRPS:2;terpene:1

    into HTML pill elements.
    """
    if (
        not isinstance(bgc_types_str, str)
        or not bgc_types_str.strip()
        or bgc_types_str.strip().upper() == "N/A"
    ):
        return "N/A"

    pills = []

    for entry in bgc_types_str.split(";"):
        entry = entry.strip()

        if not entry:
            continue

        name, separator, count = entry.rpartition(":")

        if not separator:
            name = entry
            count = "1"

        name = name.strip()
        count = count.strip()

        pills.append(
            '<span class="bgc-pill">'
            f"{safe_html(name)} &times;{safe_html(count)}"
            "</span>"
        )

    return "".join(pills) if pills else "N/A"


def load_summary(summary_path):
    """Read and normalize the BGC-per-bin summary table."""
    try:
        summary_df = pd.read_csv(summary_path, sep="\t")
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Summary file does not exist: {summary_path}"
        ) from exc
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

    if "bin_id" not in summary_df.columns:
        raise ValueError(
            f"Required column 'bin_id' is missing from {summary_path}"
        )

    text_columns = [
        "taxonomy",
        "taxonomy_lowest",
        "markermag_taxonomy",
        "markermag_taxonomy_lowest",
        "bgc_types",
    ]

    for column in text_columns:
        if column not in summary_df.columns:
            summary_df[column] = "N/A"
        else:
            summary_df[column] = summary_df[column].fillna("N/A")

    if "n_bgcs" not in summary_df.columns:
        summary_df["n_bgcs"] = 0
    else:
        summary_df["n_bgcs"] = (
            pd.to_numeric(summary_df["n_bgcs"], errors="coerce")
            .fillna(0)
            .astype(int)
        )

    return summary_df


def make_html(summary_df, sample_pool, report_prefix):
    """Generate the sortable HTML page."""
    sample_pool_html = safe_html(sample_pool)
    report_prefix = str(report_prefix).rstrip("/")

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">

    <meta
        name="viewport"
        content="width=device-width, initial-scale=1.0"
    >

    <title>
        antiSMASH Results per Bin - {sample_pool_html}
    </title>

    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
        }}

        table {{
            border-collapse: collapse;
            width: 100%;
            margin-top: 20px;
        }}

        th,
        td {{
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }}

        th {{
            background-color: #4CAF50;
            color: white;
            cursor: pointer;
            user-select: none;
            position: relative;
            padding-right: 26px;
        }}

        th:hover {{
            background-color: #449d48;
        }}

        th.sort-asc::after {{
            content: "▲";
            position: absolute;
            right: 8px;
        }}

        th.sort-desc::after {{
            content: "▼";
            position: absolute;
            right: 8px;
        }}

        tbody tr:hover {{
            background-color: #f5f5f5;
        }}

        a {{
            color: #0066cc;
            text-decoration: none;
        }}

        a:hover {{
            text-decoration: underline;
        }}

        .bgc-pill {{
            display: inline-block;
            background-color: #e8f0fe;
            color: #1a3a6b;
            border-radius: 12px;
            padding: 2px 10px;
            margin: 2px;
            font-size: 0.85em;
            white-space: nowrap;
        }}
    </style>
</head>

<body>
    <h1>antiSMASH BGC Results per MAG</h1>
    <h2>Sample Pool: {sample_pool_html}</h2>

    <table id="antismash-table">
        <thead>
            <tr>
                <th data-type="text">MAG ID</th>
                <th data-type="text">BAT Taxonomy</th>
                <th data-type="text">MarkerMAG Taxonomy</th>
                <th data-type="number">Number of BGCs</th>
                <th data-type="text">BGC Types</th>
                <th data-type="text">antiSMASH Report</th>
            </tr>
        </thead>

        <tbody>
"""

    if summary_df.empty:
        html_content += """
            <tr>
                <td colspan="6">No bins with BGCs found</td>
            </tr>
"""
    else:
        for _, row in summary_df.iterrows():
            bin_id_raw = str(row["bin_id"])
            bin_id_html = safe_html(bin_id_raw)

            taxonomy_full = safe_html(
                row.get("taxonomy", "N/A")
            )
            taxonomy_lowest = safe_html(
                row.get("taxonomy_lowest", "N/A")
            )

            markermag_full = safe_html(
                row.get("markermag_taxonomy", "N/A")
            )
            markermag_lowest = safe_html(
                row.get("markermag_taxonomy_lowest", "N/A")
            )

            n_bgcs = int(row.get("n_bgcs", 0))

            bgc_types_html = bgc_types_to_pills(
                row.get("bgc_types", "N/A")
            )

            report_href = (
                f"{report_prefix}/"
                f"{bin_id_html}/index.html"
            )

            html_content += f"""
            <tr>
                <td>{bin_id_html}</td>

                <td>
                    <em title="{taxonomy_full}">
                        {taxonomy_lowest}
                    </em>
                </td>

                <td>
                    <em title="{markermag_full}">
                        {markermag_lowest}
                    </em>
                </td>

                <td>{n_bgcs}</td>

                <td>{bgc_types_html}</td>

                <td>
                    <a
                        href="{report_href}"
                        target="_blank"
                        rel="noopener noreferrer"
                    >
                        View Report
                    </a>
                </td>
            </tr>
"""

    html_content += """
        </tbody>
    </table>

    <script>
        document.addEventListener("DOMContentLoaded", function () {
            const table =
                document.getElementById("antismash-table");

            const headers =
                table.querySelectorAll("thead th");

            const tbody =
                table.querySelector("tbody");

            headers.forEach(function (header, columnIndex) {
                header.addEventListener("click", function () {
                    const rows = Array.from(
                        tbody.querySelectorAll("tr")
                    );

                    if (
                        rows.length === 1 &&
                        rows[0].querySelector("td[colspan]")
                    ) {
                        return;
                    }

                    const dataType =
                        header.dataset.type || "text";

                    const ascending =
                        !header.classList.contains("sort-asc");

                    headers.forEach(function (otherHeader) {
                        otherHeader.classList.remove(
                            "sort-asc",
                            "sort-desc"
                        );
                    });

                    header.classList.add(
                        ascending
                            ? "sort-asc"
                            : "sort-desc"
                    );

                    rows.sort(function (rowA, rowB) {
                        const valueA =
                            rowA.cells[columnIndex]
                                .textContent
                                .trim();

                        const valueB =
                            rowB.cells[columnIndex]
                                .textContent
                                .trim();

                        const missingA =
                            valueA === "" ||
                            valueA.toUpperCase() === "N/A";

                        const missingB =
                            valueB === "" ||
                            valueB.toUpperCase() === "N/A";

                        // Keep missing values at the bottom.
                        if (missingA && missingB) {
                            return 0;
                        }

                        if (missingA) {
                            return 1;
                        }

                        if (missingB) {
                            return -1;
                        }

                        let comparison;

                        if (dataType === "number") {
                            comparison =
                                Number(valueA) - Number(valueB);
                        } else {
                            comparison = valueA.localeCompare(
                                valueB,
                                undefined,
                                {
                                    numeric: true,
                                    sensitivity: "base"
                                }
                            );
                        }

                        return ascending
                            ? comparison
                            : -comparison;
                    });

                    rows.forEach(function (row) {
                        tbody.appendChild(row);
                    });
                });
            });
        });
    </script>
</body>
</html>
"""

    return html_content


def main():
    summary_path = Path(snakemake.input.bin_taxonomy)
    output_path = Path(snakemake.output.index)

    sample_pool = str(snakemake.params.sample_pool)
    report_prefix = str(snakemake.params.report_prefix)

    summary_df = load_summary(summary_path)

    html_content = make_html(
        summary_df=summary_df,
        sample_pool=sample_pool,
        report_prefix=report_prefix,
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path.write_text(
        html_content,
        encoding="utf-8",
    )

    print(f"Generated HTML index: {output_path}")
    print(f"Sample pool: {sample_pool}")
    print(f"Report prefix: {report_prefix}")
    print(f"Number of MAG rows: {len(summary_df)}")


if __name__ == "__main__":
    main()