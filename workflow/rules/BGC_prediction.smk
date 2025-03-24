""" The rules related to Biosynthetic Gene Cluster (BGC) prediction and related analyses"""

rule antismash:
    input:
        gff=f"{outdir}/results/04_gene_prediction/augustify/{{sample_pool}}/{{sample_pool}}_eukproteins.gff"

    output:
        html=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/antismash.html",
        json=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/antismash.json"
    conda:
        "../envs/antismash.yaml"
    params:
        outdir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial",
        threads=config['threads'],
        database_dir=config['antismash_db']

    shell:
        "antismash -c {params.threads} --genefinding-gff3 {input.gff} --output-dir {params.outdir} \
        --taxon bacteria --output-basename bacterial --cc-mibig --cb-general --cb-knownclusters --databases {params.database_dir}"

rule fungismash:
    input:
        gff=f"{outdir}/results/04_gene_prediction/augustify/{{sample_pool}}/{{sample_pool}}_eukproteins.gff"

    output:
        html=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/antismash.html",
        json=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/antismash.json"
    conda:
        "../envs/antismash.yaml"
    params:
        outdir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/",
        threads=config['threads'],
        database_dir=config['antismash_db']

    shell:
        "antismash -c {params.threads} --genefinding-gff3 {input.gff} --output-dir {params.outdir} \
        --taxon fungi --cassis --output-basename fungal --cc-mibig --cb-general --cb-knownclusters --databases {params.database_dir}"

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
