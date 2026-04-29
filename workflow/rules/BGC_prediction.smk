""" The rules related to Biosynthetic Gene Cluster (BGC) prediction and related analyses"""

minsize_antismash = config['min_contig_antismash']
antismash_db = config["antismash_db"]
download_antismash_db = config.get("download_antismash_db", False)


def antismash_db_marker(wildcards):
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
    shell:
        """
        mkdir -p {params.dbdir}
        download-antismash-databases --output-dir {params.dbdir}
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

    shell:
        """
        test -d {params.database_dir} || \ 
        (echo "ERROR: antiSMASH database directory not found: {params.database_dir}") >&2; exit 1) 
        
        antismash {input.contigs} -c {params.threads} --genefinding-gff3 {input.gff} --output-dir {params.outdir} \
        --taxon bacteria --output-basename bacterial --cc-mibig --cb-general --cb-knownclusters --databases {params.database_dir}
        """

rule fungismash:
    input:
        db_ready=antismash_db_input,
        gff=f"{outdir}/results/04_gene_prediction/augustify/{{sample_pool}}/{{sample_pool}}_eukproteins.gff",
        contigs=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample_pool}}/eukaryotes.fasta"

    output:
        html=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/index.html",
        json=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/fungal.json"
    conda:
        "../envs/antismash.yaml"
    params:
        outdir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/",
        threads=config['threads'],
        database_dir=antismash_db

    shell:
        """
        test -d {params.database_dir} || \ 
        (echo "ERROR: antiSMASH database directory not found: {params.database_dir}") >&2; exit 1) 
        
        antismash {input.contigs} -c {params.threads} --genefinding-gff3 {input.gff} --output-dir {params.outdir} \
        --taxon fungi --cassis --output-basename fungal --cc-mibig --cb-general --cb-knownclusters --databases {params.database_dir}
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
