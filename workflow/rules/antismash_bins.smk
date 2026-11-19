rule map_bgc_to_bins:
    input:
        bgc_summary=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/bacterial_summary.tsv",
        contig2bin=f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool_contig2bin.tsv",
        bat_taxonomy=f"{outdir}/results/06_binning/BAT/{{sample_pool}}/{{sample_pool}}.bin2classification.names.txt",
        markermag_taxonomy=f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag/{{sample_pool}}_linkages_by_genome_taxonomy.txt"
    output:
        bgc_bin_taxonomy=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bgc_bin_taxonomy.tsv"
    conda:
        "../envs/python_simple.yaml"
    log:
        f"{outdir}/logs/map_bgc_to_bins/{{sample_pool}}.log"
    shell:
        """
        python workflow/scripts/map_bgc_to_bins.py \
            --bgc {input.bgc_summary:q} \
            --contig2bin {input.contig2bin:q} \
            --bat {input.bat_taxonomy:q} \
            --markermag {input.markermag_taxonomy:q} \
            --output {output.bgc_bin_taxonomy:q} \
            > {log:q} 2>&1
        """


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
        html_files=aggregate_bin_htmls,
        bin_taxonomy=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bgc_per_bin_summary.tsv"
    output:
        index=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/per_bin_index.html"
    params:
        sample_pool=lambda wildcards: wildcards.sample_pool,

        # Path used after per_bin_index.html is copied to 09_plots/
        report_prefix="per_bin"
    conda:
        "../envs/python_simple.yaml"
    script:
        "../scripts/make_per_bin_antismash_index.py"