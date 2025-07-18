""" Rules related to gene prediction """

rule prodigal:
    input:
        contigs=branch(config['assembly_method'] == "coassembly",
            then=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample}}/{{sample}}_assembly.fasta",
            otherwise=f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}/{{sample}}_assembly.fasta")

    output:
        gff=f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}/{{sample}}_genes.gff",
        faa=f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}/{{sample}}_proteins.faa",
        fna=f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}/{{sample}}_orfs.fna"

    conda:
        "../envs/prodigal.yaml"

    shell:"prodigal -i {input.contigs} -o {output.gff} -a {output.faa} -d {output.fna} -p meta -f gff"

rule whokaryote:
    input:
        contigs=branch(config['assembly_method'] == "coassembly",
            then=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample}}/{{sample}}_assembly.fasta",
            otherwise=f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}/{{sample}}_assembly.fasta"),
        prodigal_gff=f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}/{{sample}}_genes.gff"

    output:
        headers_euk=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/eukaryote_contig_headers.txt",
        headers_prok=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/prokaryote_contig_headers.txt",
        contigs_size=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/contigs1000.fasta",
        euk_fasta=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/eukaryotes.fasta",
        prok_fasta=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/prokaryotes.fasta",
        unclassified_fasta=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/unclassified.fasta"

    conda:
        "../envs/whokaryote.yaml"

    params:
        outdir=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}"

    shell:
        "whokaryote.py --contigs {input.contigs} --outdir {params.outdir} --prodigal_file {input.prodigal_gff} --minsize 1000 --f"


rule augustify:
    input:
        f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/eukaryotes.fasta"
    output:
        genes=f"{outdir}/results/04_gene_prediction/augustify/{{sample}}/{{sample}}_eukgenes.gff",
        proteins=f"{outdir}/results/04_gene_prediction/augustify/{{sample}}/{{sample}}_eukproteins.faa",
    conda:
        "../envs/augustus.yaml"
    params:
        param_file="../resources/augustify_params.txt",
        threads=config['threads']
    script:
        "../scripts/augustify.py -g {input} -p {params.param_file} -m {output.genes} -P {output.proteins} -t {params.threads}"


