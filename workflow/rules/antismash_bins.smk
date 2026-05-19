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

        # Load BAT taxonomy (skip comment lines starting with #)
        bat_df = pd.read_csv(input.bat_taxonomy,sep="\t",comment="#")
        # BAT output has columns: bin, classification, reason, lineage, lineage scores
        bat_df = bat_df.rename(columns={bat_df.columns[0]: "bin_id"})

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

rule summarize_bgc_per_bin:
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

def get_bins_for_sample_pool(wildcards):
    """Get list of bins for a sample pool from the contig2bin file."""
    import pandas as pd
    
    contig2bin_file = checkpoints.dastool.get(
        sample_pool=wildcards.sample_pool
    ).output.c2bin
    
    # Read contig2bin file
    df = pd.read_csv(contig2bin_file, sep="\t", header=None, names=["contig", "bin"])
    bins = df["bin"].unique().tolist()
    
    return bins


# checkpoint list_bins:
#     """Create a file listing all bins for a sample pool."""
#     input:
#         contig2bin=f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool_contig2bin.tsv"
#     output:
#         bin_list=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/bin_list.txt"
#     run:
#         import pandas as pd
#         
#         df = pd.read_csv(input.contig2bin, sep="\t", header=None, names=["contig", "bin"])
#         bins = sorted(df["bin"].unique())
#         
#         with open(output.bin_list, "w") as f:
#             for bin_id in bins:
#                 f.write(f"{bin_id}\n")


def aggregate_bin_htmls(wildcards):
    """Aggregate all bin HTML files for a sample pool."""
    import pandas as pd
    
    # Get the contig2bin file from the dastool checkpoint
    checkpoint_output = checkpoints.dastool.get(
        sample_pool=wildcards.sample_pool
    ).output.c2bin
    
    # Read the contig2bin file to get list of bins
    try:
        df = pd.read_csv(checkpoint_output, sep="\t", header=None, names=["contig", "bin"])
        bins = sorted(df["bin"].unique().tolist())
    except (FileNotFoundError, pd.errors.EmptyDataError):
        # During dry-run or if file is empty
        bins = []
    
    # Return list of expected HTML files
    return expand(
        f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin/{{bin_id}}/index.html",
        sample_pool=wildcards.sample_pool,
        bin_id=bins
    )


rule filter_antismash_by_bin:
    """Filter antiSMASH GenBank files to only include contigs from one bin."""
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


rule regenerate_antismash_html_per_bin:
    """Regenerate antiSMASH HTML from filtered GenBank files for a bin."""
    input:
        filtered_gbk=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin/{{bin_id}}_filtered_gbk"
    output:
        html=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin/{{bin_id}}/index.html"
    params:
        output_dir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin/{{bin_id}}",
        bin_id="{bin_id}"
    threads: 1
    conda:
        "../envs/antismash.yaml"
    script:
        "../scripts/regenerate_antismash_html.py"

# rule regenerate_antismash_html_per_bin:
#     """Regenerate antiSMASH HTML from filtered GenBank files for a bin."""
#     input:
#         filtered_gbk=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin/{{bin_id}}_filtered_gbk"
#     output:
#         html=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin/{{bin_id}}/index.html"
#     params:
#         output_dir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin/{{bin_id}}",
#         bin_id="{bin_id}",
#         database_dir=antismash_db
#     threads: 1
#     conda:
#         "../envs/antismash.yaml"
#     shell:
#         """
#         # Check if there are any GenBank files
#         n_gbk=$(find {input.filtered_gbk} -name "*.gbk" | wc -l)
#
#         if [ "$n_gbk" -eq 0 ]; then
#             # Create empty HTML if no BGCs in this bin
#             mkdir -p {params.output_dir}
#             echo "<html><body><h1>No BGCs found in bin {params.bin_id}</h1></body></html>" > {output.html}
#         else
#             # Copy GenBank files to output directory
#             mkdir -p {params.output_dir}
#             cp {input.filtered_gbk}/*.gbk {params.output_dir}/
#
#             # Use antiSMASH to regenerate HTML from GenBank files
#             # This is much faster than re-running full analysis
#             cd {params.output_dir}
#
#             # Create a minimal completion marker so antiSMASH knows these are pre-analyzed
#             # Then just generate the HTML visualization
#             for gbk in *.gbk; do
#                 python -c "
# from Bio import SeqIO
# import json
#
# # antiSMASH can regenerate HTML if the GenBank files are properly formatted
# # The files from the original run should already have all the necessary features
# print('GenBank file ready: $gbk')
#                 "
#             done
#         fi
#         """


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


