rule map_bgc_to_bins:
    input:
        bgc_summary=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/bacterial_summary.tsv",
        contig2bin=f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool_contig2bin.tsv",
        bat_taxonomy=f"{outdir}/results/06_binning/BAT/{{sample_pool}}/{{sample_pool}}.bin2classification.names.txt",
        markermag_taxonomy=f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag/{{sample_pool}}_linkages_by_genome_taxonomy.txt"
    output:
        bgc_bin_taxonomy=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bgc_bin_taxonomy.tsv"
    run:
        import pandas as pd

        # Load BGC summary
        bgc_df = pd.read_csv(input.bgc_summary,sep="\t")

        # Load contig to bin mapping
        c2b_df = pd.read_csv(input.contig2bin,sep="\t",header=None,names=["contig_id", "bin_id"])

        # Load BAT taxonomy. The header line itself starts with "#" (e.g.
        # "# bin\tclassification\treason\t..."), so comment="#" must NOT be
        # used here -- it would drop the header and treat the first data row
        # as column names instead.
        bat_df = pd.read_csv(input.bat_taxonomy,sep="\t")
        bat_df = bat_df.rename(columns={bat_df.columns[0]: "bin_id"})
        # CAT_pack bins was run with -s .fa, but add_names does not strip
        # that suffix from the bin name (e.g. "A_bin.001.fa"), while
        # DASTool_contig2bin.tsv uses the bare bin id (e.g. "A_bin.001").
        # Strip it so the merge below on bin_id actually matches.
        bat_df["bin_id"] = bat_df["bin_id"].str.removesuffix(".fa")

        # CAT_pack's own "lineage" column is a taxid string (e.g.
        # "1;131567;2;..."), not a human-readable name. Build a readable
        # lineage from the named rank columns added by --only_official,
        # skipping ranks with no support/NA.
        rank_cols = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
        bat_df["lineage"] = bat_df[rank_cols].apply(
            lambda row: ";".join(
                str(v) for v in row if pd.notna(v) and v not in ("no support", "NA")
            ),
            axis=1
        )
        # Lowest available rank, kept alongside the full lineage -- the full
        # string stays in the tsv for downstream figure scripts, while the
        # HTML report only has room to show the lowest rank per bin.
        bat_df["lineage_lowest"] = bat_df["lineage"].apply(
            lambda s: s.split(";")[-1] if s else "N/A"
        )

        # Merge BGC with bin assignment
        bgc_bins = bgc_df.merge(c2b_df,on="contig_id",how="left")

        # Merge with taxonomy
        bgc_bins_tax = bgc_bins.merge(
            bat_df[["bin_id", "lineage", "lineage_lowest"]],
            on="bin_id",
            how="left"
        )

        # Load markerMAG taxonomy (16S-based genome-marker linkage). Not every
        # bin is reached by a 16S linkage, so this is a left join and missing
        # values are expected -- shown alongside, not instead of, BAT taxonomy.
        try:
            markermag_df = pd.read_csv(input.markermag_taxonomy, sep="\t")
        except (FileNotFoundError, pd.errors.EmptyDataError):
            markermag_df = pd.DataFrame(columns=["GenomicSeq", "taxonomy"])

        if not markermag_df.empty and "GenomicSeq" in markermag_df.columns:
            markermag_df = markermag_df.rename(
                columns={"GenomicSeq": "bin_id", "taxonomy": "markermag_taxonomy"}
            )
            # "unknown" is add_markermag_taxonomy.py's placeholder for marker
            # genes with no phyloFlash classification -- not informative here.
            markermag_df = markermag_df[markermag_df["markermag_taxonomy"] != "unknown"]
            markermag_df["markermag_taxonomy_lowest"] = markermag_df["markermag_taxonomy"].apply(
                lambda s: s.split(";")[-1] if isinstance(s, str) and s else "N/A"
            )
            # A bin can have multiple marker gene linkages (e.g. several 16S
            # copies); collapse to one row per bin. "||" separates distinct
            # full lineages (each internally ";"-delimited by rank) so it
            # stays unambiguous; the lowest-rank names use ", " since they're
            # single tokens.
            markermag_tax = markermag_df.groupby("bin_id").agg(
                markermag_taxonomy=("markermag_taxonomy", lambda x: " || ".join(sorted(set(x)))),
                markermag_taxonomy_lowest=("markermag_taxonomy_lowest", lambda x: ", ".join(sorted(set(x))))
            ).reset_index()
        else:
            markermag_tax = pd.DataFrame(columns=["bin_id", "markermag_taxonomy", "markermag_taxonomy_lowest"])

        bgc_bins_tax = bgc_bins_tax.merge(markermag_tax, on="bin_id", how="left")

        # Add indicator for unbinned contigs
        bgc_bins_tax["bin_status"] = bgc_bins_tax["bin_id"].apply(
            lambda x: "binned" if pd.notna(x) else "unbinned"
        )

        # Save combined table
        bgc_bins_tax.to_csv(output.bgc_bin_taxonomy,sep="\t",index=False)


rule combine_all_bgc_bin_taxonomy:
    input:
        bgc_bin_tax_files=expand(
            f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bgc_bin_taxonomy.tsv",
            sample_pool=ASSEMBLY_UNITS
        )
    output:
        combined=f"{outdir}/results/08_BGC/antismash/combined_bgc_bin_taxonomy.tsv"
    run:
        import pandas as pd

        dfs = []
        for file in input.bgc_bin_tax_files:
            df = pd.read_csv(file,sep="\t")
            dfs.append(df)

        combined_df = pd.concat(dfs,ignore_index=True)
        combined_df.to_csv(output.combined,sep="\t",index=False)

checkpoint summarize_bgc_per_bin:
    """Summarize BGCs per bin. A checkpoint because downstream per-bin
    antiSMASH report generation only runs for bins listed here, i.e. bins
    with at least one detected BGC -- the exact bin list isn't known until
    this rule has run."""
    input:
        bgc_bin_tax=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bgc_bin_taxonomy.tsv"
    output:
        bin_summary=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bgc_per_bin_summary.tsv"
    run:
        import pandas as pd

        df = pd.read_csv(input.bgc_bin_tax,sep="\t")

        # Only binned BGCs
        binned = df[df["bin_status"] == "binned"].copy()

        summary_columns = [
            "bin_id", "taxonomy", "taxonomy_lowest",
            "markermag_taxonomy", "markermag_taxonomy_lowest",
            "n_bgcs", "bgc_types", "sample_pool"
        ]

        if len(binned) == 0:
            # Create empty output with expected columns
            summary = pd.DataFrame(columns=summary_columns)
        else:
            # Group by bin and summarize. Taxonomy (both full and
            # lowest-rank) is the same for all BGCs in a bin.
            summary = binned.groupby("bin_id").agg({
                "lineage": "first",
                "lineage_lowest": "first",
                "markermag_taxonomy": "first",
                "markermag_taxonomy_lowest": "first",
                "bgc_id": "count",  # Count BGCs
                # Unique BGC types with their per-bin counts, e.g.
                # "NRPS:2;T1PKS:1", most frequent first -- lets the HTML
                # report render these as pills without losing the counts.
                "bgc_product": lambda x: ";".join(
                    f"{name}:{count}" for name, count in
                    sorted(x.value_counts().items(), key=lambda kv: (-kv[1], kv[0]))
                ),
                "sample_pool": "first"
            }).reset_index()

            summary.columns = summary_columns

        summary.to_csv(output.bin_summary,sep="\t",index=False)

def bins_with_bgcs(wildcards):
    """Bin IDs that have at least one detected BGC for this sample pool.

    Reads the summarize_bgc_per_bin checkpoint, which only lists bins with
    binned BGC regions -- bins with zero BGCs simply don't appear here, so
    no filtering/regeneration work is done for them.
    """
    import pandas as pd

    summary_file = checkpoints.summarize_bgc_per_bin.get(
        sample_pool=wildcards.sample_pool
    ).output.bin_summary

    try:
        df = pd.read_csv(summary_file, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return []

    if df.empty or "bin_id" not in df.columns:
        return []

    return sorted(df["bin_id"].unique().tolist())


def aggregate_bin_htmls(wildcards):
    """Expected per-bin antiSMASH HTML report, for every bin with a BGC."""
    bins = bins_with_bgcs(wildcards)
    return expand(
        f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin/{{bin_id}}/index.html",
        sample_pool=wildcards.sample_pool,
        bin_id=bins
    )


rule filter_antismash_by_bin:
    """Filter antiSMASH GenBank region files to only include contigs from one
    bin. Not used by the per-bin HTML report (see filter_antismash_json_by_bin
    for that) -- kept as an input source for tools that consume per-bin
    GenBank files directly, e.g. BiG-SCAPE."""
    input:
        antismash_json=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/bacterial.json",
        contig2bin=f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool_contig2bin.tsv"
    output:
        filtered_dir=directory(f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin/{{bin_id}}_filtered_gbk")
    params:
        antismash_dir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial",
        bin_id="{bin_id}"
    conda:
        "../envs/antismash.yaml"
    script:
        "../scripts/filter_antismash_bins.py"


rule filter_antismash_json_by_bin:
    """Filter the antiSMASH results JSON down to the records belonging to
    one bin. This is the input antiSMASH's own --reuse-results mode needs to
    regenerate a per-bin HTML report without re-running any analysis."""
    input:
        antismash_json=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/bacterial.json",
        contig2bin=f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool_contig2bin.tsv"
    output:
        filtered_json=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin/{{bin_id}}_reuse.json"
    params:
        bin_id="{bin_id}"
    log:
        f"{outdir}/logs/antismash/bacterial/per_bin/{{sample_pool}}_{{bin_id}}_filter_json.log"
    script:
        "../scripts/filter_antismash_json_by_bin.py"


rule regenerate_antismash_html_per_bin:
    """Regenerate a full antiSMASH HTML report for one bin using antiSMASH's
    own --reuse-results mode against the bin's filtered results JSON. Only
    output rendering is redone; gene finding, cluster detection, ClusterBlast
    etc. are all skipped since their results are already cached in the JSON,
    so this is much cheaper than a full antiSMASH run."""
    input:
        filtered_json=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin/{{bin_id}}_reuse.json"
    output:
        html=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin/{{bin_id}}/index.html"
    params:
        output_dir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin/{{bin_id}}",
        bin_id="{bin_id}",
        database_dir=antismash_db
    threads: 2
    log:
        f"{outdir}/logs/antismash/bacterial/per_bin/{{sample_pool}}_{{bin_id}}_regenerate_html.log"
    conda:
        "../envs/antismash.yaml"
    shell:
        """
        antismash --reuse-results {input.filtered_json} \
            --output-dir {params.output_dir} \
            --output-basename {params.bin_id} \
            --databases {params.database_dir} \
            -c {threads} \
            > {log} 2>&1
        """


rule aggregate_bin_antismash_reports:
    """Create a sortable index page linking to all per-bin antiSMASH reports."""
    input:
        bin_taxonomy=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bgc_per_bin_summary.tsv"
    output:
        index=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin_index.html"
    params:
        report_prefix="per_bin",
        sample_pool="{sample_pool}"
    run:
        import pandas as pd
        import html
        from pathlib import Path
        # Load the per-bin summary table.
        try:
            summary_df = pd.read_csv(
                input.bin_taxonomy,
                sep="\t"
            )
        except (FileNotFoundError, pd.errors.EmptyDataError):
            summary_df = pd.DataFrame()

        def safe_html(value):
            """Escape a value before inserting it into HTML."""
            if pd.isna(value):
                return "N/A"

            return html.escape(str(value), quote=True)

        def bgc_types_to_pills(bgc_types_str):
            """
            Convert:
                NRPS:2;terpene:1

            into styled HTML pill elements.
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

        if not summary_df.empty:
            required_columns = ["bin_id"]

            missing_columns = [
                column
                for column in required_columns
                if column not in summary_df.columns
            ]

            if missing_columns:
                raise ValueError(
                    "Missing required column(s) in "
                    f"{input.bin_taxonomy}: "
                    + ", ".join(missing_columns)
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
                    summary_df[column] = (
                        summary_df[column]
                        .fillna("N/A")
                        .astype(str)
                    )

            if "n_bgcs" not in summary_df.columns:
                summary_df["n_bgcs"] = 0
            else:
                summary_df["n_bgcs"] = (
                    pd.to_numeric(
                        summary_df["n_bgcs"],
                        errors="coerce"
                    )
                    .fillna(0)
                    .astype(int)
                )

        sample_pool_html = safe_html(params.sample_pool)

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
            color: #222;
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
            vertical-align: middle;
        }}

        th {{
            background-color: #4CAF50;
            color: white;
            cursor: pointer;
            user-select: none;
            position: relative;
            padding-right: 28px;
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

    <h2>
        Sample Pool: {sample_pool_html}
    </h2>

    <table id="antismash-table">
        <thead>
            <tr>
                <th data-type="text">
                    MAG ID
                </th>

                <th data-type="text">
                    BAT Taxonomy
                </th>

                <th data-type="text">
                    MarkerMAG Taxonomy
                </th>

                <th data-type="number">
                    Number of BGCs
                </th>

                <th data-type="text">
                    BGC Types
                </th>

                <th data-type="text">
                    antiSMASH Report
                </th>
            </tr>
        </thead>

        <tbody>
"""

        if summary_df.empty:
            html_content += """
            <tr>
                <td colspan="6">
                    No bins with BGCs found
                </td>
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
                    row.get(
                        "markermag_taxonomy_lowest",
                        "N/A"
                    )
                )

                n_bgcs = int(row.get("n_bgcs", 0))

                bgc_types_html = bgc_types_to_pills(
                    row.get("bgc_types", "N/A")
                )

                report_href = (
                    f"{params.report_prefix}/"
                    f"{bin_id_html}/index.html"
                )

                html_content += f"""
            <tr>
                <td>
                    {bin_id_html}
                </td>

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

                <td>
                    {n_bgcs}
                </td>

                <td>
                    {bgc_types_html}
                </td>

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
        document.addEventListener(
            "DOMContentLoaded",
            function () {
                const table = document.getElementById(
                    "antismash-table"
                );

                const headers = table.querySelectorAll(
                    "thead th"
                );

                const tbody = table.querySelector(
                    "tbody"
                );

                headers.forEach(
                    function (header, columnIndex) {
                        header.addEventListener(
                            "click",
                            function () {
                                const rows = Array.from(
                                    tbody.querySelectorAll(
                                        "tr"
                                    )
                                );

                                /*
                                 * Do not sort the single
                                 * "no results" row.
                                 */
                                if (
                                    rows.length === 1 &&
                                    rows[0].querySelector(
                                        "td[colspan]"
                                    )
                                ) {
                                    return;
                                }

                                const dataType =
                                    header.dataset.type ||
                                    "text";

                                /*
                                 * First click: ascending.
                                 * Second click: descending.
                                 */
                                const ascending =
                                    !header.classList.contains(
                                        "sort-asc"
                                    );

                                headers.forEach(
                                    function (otherHeader) {
                                        otherHeader.classList
                                            .remove(
                                                "sort-asc",
                                                "sort-desc"
                                            );
                                    }
                                );

                                header.classList.add(
                                    ascending
                                        ? "sort-asc"
                                        : "sort-desc"
                                );

                                rows.sort(
                                    function (rowA, rowB) {
                                        const cellA =
                                            rowA.cells[
                                                columnIndex
                                            ];

                                        const cellB =
                                            rowB.cells[
                                                columnIndex
                                            ];

                                        const valueA = cellA
                                            ? cellA.textContent
                                                .trim()
                                            : "";

                                        const valueB = cellB
                                            ? cellB.textContent
                                                .trim()
                                            : "";

                                        const missingA =
                                            valueA === "" ||
                                            valueA
                                                .toUpperCase()
                                                === "N/A";

                                        const missingB =
                                            valueB === "" ||
                                            valueB
                                                .toUpperCase()
                                                === "N/A";

                                        /*
                                         * Missing values are
                                         * always placed last.
                                         */
                                        if (
                                            missingA &&
                                            missingB
                                        ) {
                                            return 0;
                                        }

                                        if (missingA) {
                                            return 1;
                                        }

                                        if (missingB) {
                                            return -1;
                                        }

                                        let comparison = 0;

                                        if (
                                            dataType ===
                                            "number"
                                        ) {
                                            comparison =
                                                Number(valueA) -
                                                Number(valueB);
                                        } else {
                                            comparison =
                                                valueA
                                                    .localeCompare(
                                                        valueB,
                                                        undefined,
                                                        {
                                                            numeric:
                                                                true,
                                                            sensitivity:
                                                                "base"
                                                        }
                                                    );
                                        }

                                        return ascending
                                            ? comparison
                                            : -comparison;
                                    }
                                );

                                rows.forEach(
                                    function (row) {
                                        tbody.appendChild(
                                            row
                                        );
                                    }
                                );
                            }
                        );
                    }
                );
            }
        );
    </script>
</body>
</html>
"""

        output_path = Path(output.index)
        output_path.parent.mkdir(
            parents=True,
            exist_ok=True
        )

        output_path.write_text(
            html_content,
            encoding="utf-8"
        )

        print(
            f"Generated sortable antiSMASH index: "
            f"{output_path}"
        )
        print(
            f"Sample pool: {params.sample_pool}"
        )
        print(
            f"Report prefix: {params.report_prefix}"
        )
        print(
            f"Number of MAGs: {len(summary_df)}"
        )