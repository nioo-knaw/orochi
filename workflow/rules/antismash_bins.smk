rule map_bgc_to_bins:
    input:
        bgc_summary=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/bacterial_summary.tsv",
        contig2bin=f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool_contig2bin.tsv",
        bat_taxonomy=f"{outdir}/results/06_binning/BAT/{{sample_pool}}/{{sample_pool}}.bin2classification.names.txt"
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

        # Merge BGC with bin assignment
        bgc_bins = bgc_df.merge(c2b_df,on="contig_id",how="left")

        # Merge with taxonomy
        bgc_bins_tax = bgc_bins.merge(
            bat_df[["bin_id", "lineage"]],
            on="bin_id",
            how="left"
        )

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

        if len(binned) == 0:
            # Create empty output with expected columns
            summary = pd.DataFrame(columns=[
                "bin_id", "taxonomy", "n_bgcs", "bgc_types", "sample_pool"
            ])
        else:
            # Group by bin and summarize
            summary = binned.groupby("bin_id").agg({
                "lineage": "first",  # Taxonomy is the same for all BGCs in a bin
                "bgc_id": "count",  # Count BGCs
                "bgc_product": lambda x: ";".join(sorted(set(x))),  # Unique BGC types
                "sample_pool": "first"
            }).reset_index()

            summary.columns = ["bin_id", "taxonomy", "n_bgcs", "bgc_types", "sample_pool"]

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
    """Create an index page linking to all per-bin antiSMASH reports."""
    input:
        html_files=aggregate_bin_htmls,
        bin_taxonomy=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bgc_per_bin_summary.tsv"
    output:
        index=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin_index.html"
    params:
        sample_pool="{sample_pool}"
    run:
        import pandas as pd
        from pathlib import Path
        
        # Load bin summary
        try:
            summary_df = pd.read_csv(input.bin_taxonomy, sep="\t")
        except (FileNotFoundError, pd.errors.EmptyDataError):
            summary_df = pd.DataFrame()
        
        # Create HTML index
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>antiSMASH Results per Bin - {params.sample_pool}</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                table {{ border-collapse: collapse; width: 100%; margin-top: 20px; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #4CAF50; color: white; }}
                tr:hover {{ background-color: #f5f5f5; }}
                a {{ color: #0066cc; text-decoration: none; }}
                a:hover {{ text-decoration: underline; }}
            </style>
        </head>
        <body>
            <h1>antiSMASH BGC Results per Bin</h1>
            <h2>Sample Pool: {params.sample_pool}</h2>
            <table>
                <tr>
                    <th>Bin ID</th>
                    <th>Taxonomy</th>
                    <th>Number of BGCs</th>
                    <th>BGC Types</th>
                    <th>antiSMASH Report</th>
                </tr>
        """
        
        if not summary_df.empty:
            for _, row in summary_df.iterrows():
                bin_id = row['bin_id']
                taxonomy = row.get('taxonomy', 'N/A')
                n_bgcs = row.get('n_bgcs', 0)
                bgc_types = row.get('bgc_types', 'N/A')
                
                html_content += f"""
                <tr>
                    <td>{bin_id}</td>
                    <td><em>{taxonomy}</em></td>
                    <td>{n_bgcs}</td>
                    <td>{bgc_types}</td>
                    <td><a href="per_bin/{bin_id}/index.html" target="_blank">View Report</a></td>
                </tr>
                """
        else:
            html_content += "<tr><td colspan='5'>No bins with BGCs found</td></tr>"
        
        html_content += """
            </table>
        </body>
        </html>
        """
        
        with open(output.index, 'w') as f:
            f.write(html_content)


