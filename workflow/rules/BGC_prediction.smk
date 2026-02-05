""" The rules related to Biosynthetic Gene Cluster (BGC) prediction and related analyses"""

rule antismash:
    # Sometime, MAG genomes can't matach the whole gff annotations. Run prodigal-m again in antismash for prokaryotic genoomes
    input:
        derep_ok = ancient(f"{outdir}/results/06_binning/drep/dereplicated_genomes/drep.done"),
        contigs = lambda wc: f"{outdir}/results/06_binning/drep/dereplicated_genomes/{wc.genome}.fa",
        # gff = f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_prokaryote_{minsize}.gff"
    output:
        done = touch(f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/{{genome}}/.antismash.done")
    params:
        outdir = f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/{{genome}}",
        threads = config["threads"],
        database_dir = config["antismash_db"],
    threads: config["threads"]
    conda:
        "../envs/antismash.yaml"
    shell:
        """
        rm -rf {params.outdir}
        antismash {input.contigs} \
          -c {threads} \
          --output-dir {params.outdir} \
          --taxon bacteria \
          --cc-mibig \
          --cb-knownclusters \
          --databases {params.database_dir} \
          --asf \
          --clusterhmmer \
          --tfbs \
          --rre \
          --allow-long-headers \
          --no-zip-output \
          --genefinding-tool prodigal-m
        """

rule antismash_all:
    input:
        lambda wc: expand(
            f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/{{genome}}/.antismash.done",
            sample_pool = wc.sample_pool,
            genome = get_drep_genomes(wc)
        )

rule fungismash:
    input:
        gff=f"{outdir}/results/04_gene_prediction/augustify/{{sample_pool}}/{{sample_pool}}_eukproteins_{minsize}.gff",
        contigs=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample_pool}}/eukaryotes.fasta"

    output:
        html=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/index.html",
        json=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/fungal.json"
    conda:
        "../envs/antismash.yaml"
    params:
        outdir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/",
        threads=config['threads'],
        database_dir=config['antismash_db']

    shell:
        """
        rm -rf {params.outdir}
        antismash {input.contigs} -c {params.threads} --genefinding-gff3 {input.gff} --output-dir {params.outdir} \
        --taxon fungi --cassis --output-basename fungal --cc-mibig --cb-knownclusters --databases {params.database_dir} \
        --asf --clusterhmmer --tfbs --rre --allow-long-headers --no-zip-output --genefinding-tool none
        """

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
