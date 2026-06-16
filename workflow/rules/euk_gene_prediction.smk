""" Rules related to gene prediction """


minsize_antismash = config['min_contig_antismash']

checkpoint whokaryote:
    input:
        contigs=f"{outdir}/results/03_assembly/size_filtered/{{sample}}_{minsize}/contigs_{{sample}}_{minsize}.fasta",
        prodigal_gff=f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}/{{sample}}_genes.gff"

    output:
        headers_euk=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/eukaryote_contig_headers.txt",
        headers_prok=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/prokaryote_contig_headers.txt",
        contigs_size=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/contigs{minsize_antismash}.fasta",
        euk_fasta=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/eukaryotes.fasta",
        prok_fasta=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/prokaryotes.fasta",
        unclassified_fasta=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/unclassified.fasta"

    conda:
        "../envs/whokaryote.yaml"

    params:
        outdir=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}",
        minsize_a=config['min_contig_antismash']

    log:
        f"{outdir}/logs/whokaryote/whokaryote_{{sample}}.log"

    shell:
        """
        whokaryote.py \
            --contigs {input.contigs} \
            --outdir {params.outdir} \
            --prodigal_file {input.prodigal_gff} \
            --minsize {params.minsize_a} \
            --f \
            > {log} 2>&1

        # Make sure all declared outputs exist.
        # Some samples may have no eukaryotic or unclassified contigs.
        touch {output.headers_euk}
        touch {output.headers_prok}
        touch {output.euk_fasta}
        touch {output.prok_fasta}
        touch {output.unclassified_fasta}
        touch {output.contigs_size}
        """


checkpoint augustify:
    input:
        f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/eukaryotes.fasta"
    output:
        tax=f"{outdir}/results/04_gene_prediction/augustify/{{sample}}/{{sample}}_eukclass.txt",
        gff=f"{outdir}/results/04_gene_prediction/augustify/{{sample}}/{{sample}}_eukproteins.gff",
    conda:
        "../envs/augustus.yaml"
    params:
        param_file=os.path.abspath("resources/augustify_params.txt"),
        script=os.path.abspath("workflow/scripts/augustify.py"),
        outdir=f"{outdir}/logs/augustify/{{sample}}",
    threads:
        config['threads']
    log:
        f"{outdir}/logs/augustify/augustify_{{sample}}.log"
    shell:
        """
        python {params.script} \
            -g {input} \
            -p {params.param_file} \
            -m {output.tax} \
            -P {output.gff} \
            -t {threads} \
            --outdir {params.outdir} \
            > {log} 2>&1
        # Ensure declared outputs exist even if augustify produced no annotation
        touch {output.tax}
        touch {output.gff}
        """


# This rule filters the prodigal gff file so that only prokaryotic genes on contigs above the size threshold are kept
# This step takes place after whokaryote, and is meant to be used as antismash input
rule filter_prokaryote_gff:
    input:
        gff=f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}/{{sample}}_genes.gff",
        headers_prok=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/prokaryote_contig_headers.txt"
    output:
        f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}/{{sample}}_prokaryote_{minsize_antismash}.gff"
    conda:
        "../envs/size_filter.yaml"
    params:
        size=config['min_contig_antismash'],
        outdir=f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}",
        script=os.path.abspath("workflow/scripts/filter_annotations.py")
    log:
        f"{outdir}/logs/filter_prokaryote_gff/{{sample}}_{{minsize_antismash}}.log"
    shell:
        "python {params.script} \
            --gff {input.gff} \
            --headerfile {input.headers_prok} \
            --outdir {params.outdir} \
            --minsize {params.size} \
            --sample_name {wildcards.sample} \
            2> {log}"
