""" The rules related to Biosynthetic Gene Cluster (BGC) prediction and related analyses"""
from pathlib import Path
import json

def fungismash_input(wildcards):
    whokaryote_ckpt = checkpoints.whokaryote.get(sample=wildcards.sample_pool)
    euk_fasta = whokaryote_ckpt.output.euk_fasta

    if not fasta_has_records(euk_fasta):
        return {
            "contigs": euk_fasta,
        }

    augustify_ckpt = checkpoints.augustify.get(sample=wildcards.sample_pool)
    gff = augustify_ckpt.output.gff

    return {
        "db_ready": antismash_db_input(wildcards),
        "gff": gff,
        "contigs": euk_fasta,
    }

minsize_antismash = config['min_contig_antismash']
antismash_db = config["antismash_db"]
download_antismash_db = config.get("download_antismash_db", False)


def antismash_db_input(wildcards):
    if download_antismash_db:
        return os.path.join(antismash_db, ".download_complete")
    return []

rule download_antismash_databases:
    output:
        touch(os.path.join(antismash_db, ".download_complete"))
    conda:
        "../envs/antismash.yaml"
    params:
        dbdir=antismash_db,
    log:
        f"{outdir}/logs/antismash/download_antismash_databases.log"
    shell:
        """
        mkdir -p {params.dbdir}
        download-antismash-databases --database-dir {params.dbdir} > {log} 2>&1
        """


rule antismash:
    input:
        db_ready=antismash_db_input,
        gff=f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_prokaryote_{minsize_antismash}.gff",
        contigs=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample_pool}}/prokaryotes.fasta"
    output:
        html=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/index.html",
        json=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/bacterial.json"
    conda:
        "../envs/antismash.yaml"
    params:
        outdir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial",
        threads=config['threads'],
        database_dir=antismash_db

    log:
        f"{outdir}/logs/antismash/bacterial/antismash_{{sample_pool}}.log"

    shell:
        """
        test -d {params.database_dir:q} || \
            (echo "ERROR: antiSMASH database directory not found: {params.database_dir:q}" >&2; exit 1)

        # Safety check before removing antiSMASH output directory
        test "$(basename "{params.outdir:q}")" = "bacterial" || \
            (echo "ERROR: refusing to remove unexpected output directory: {params.outdir:q}" >&2; exit 1)

        # Clean antiSMASH folder for re-run
        rm -rf {params.outdir:q}

        antismash {input.contigs:q} \
        -c {params.threads} \
        --genefinding-gff3 {input.gff:q} \
        --output-dir {params.outdir:q} \
        --taxon bacteria \
        --output-basename bacterial \
        --tfbs \
        --cc-mibig \
        --cb-general \
        --cb-knownclusters \
        --databases {params.database_dir:q} \
        > {log} 2>&1
        """

rule fungismash:
    input:
        unpack(fungismash_input)

    output:
        html=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/index.html",
        json=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/fungal.json"

    log:
        f"{outdir}/logs/antismash/fungal/fungismash_{{sample_pool}}.log"

    conda:
        "../envs/antismash.yaml"

    params:
        outdir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/",
        threads=config['threads'],
        database_dir=antismash_db

    script:
        "../scripts/run_fungismash.py"


rule summarize_antismash_bacterial:
    input:
        json=rules.antismash.output.json
    output:
        tsv=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/bacterial_summary.tsv"
    conda:
        "../envs/antismash.yaml"
    params:
        antismash_dir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial",
        taxon="bacteria"
    log:
        f"{outdir}/logs/antismash/bacterial/summarize_antismash_bacterial/{{sample_pool}}.log"
    script:
        "../scripts/summarize_antismash.py"


rule summarize_antismash_fungal:
    input:
        json=rules.fungismash.output.json
    output:
        tsv=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/fungal_summary.tsv"
    conda:
        "../envs/antismash.yaml"
    params:
        antismash_dir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal",
        taxon="fungi"
    log:
        f"{outdir}/logs/antismash/fungal/summarize_antismash_fungal/{{sample_pool}}.log"
    script:
        "../scripts/summarize_antismash.py"


rule combine_bgc_summaries:
    """Merge the bacterial and fungal antiSMASH summaries into one table.

    One row per BGC region, annotated with CAT contig-level taxonomy (the
    only taxonomy source that covers fungal contigs too) and, where the
    contig was binned, its DASTool bin."""
    input:
        bacterial_summary=rules.summarize_antismash_bacterial.output.tsv,
        # Empty whenever the assembly unit has no eukaryotic contigs with
        # predicted genes -- fungiSMASH is not run in that case.
        fungal_summary=fungal_summary_targets,
        contig2bin=f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool_contig2bin.tsv",
        cat_taxonomy=rules.CAT.output.names
    output:
        combined_table=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/combined_bgc_table.tsv"
    params:
        sample_pool="{sample_pool}"
    conda:
        "../envs/python_simple.yaml"
    script:
        "../scripts/build_combined_bgc_table.py"


rule combined_bgc_overview:
    """Standalone HTML overview of all BGCs (bacterial + fungal) in one
    assembly unit: summary statistics, charts per type, and a sortable and
    filterable table."""
    input:
        combined_table=rules.combine_bgc_summaries.output.combined_table
    output:
        index=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/combined_bgc_index.html"
    params:
        sample_pool="{sample_pool}"
    conda:
        "../envs/python_simple.yaml"
    script:
        "../scripts/make_combined_antismash_index.py"


rule bigscape:
    input:
        "path/to/antismash_output"

    output:
        "path/to/output"

    shell:
        "bigscape -options"

rule itol_bgc:
    input:
        "path/to/input"
    output:
        "path/to/output"
