""" Rules related to gene prediction """

rule prodigal:
    input:
        contigs=f"{outdir}/results/03_assembly/size_filtered/{{sample}}_{minsize}/contigs_{{sample}}_{minsize}.fasta"

    output:
        gff=f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}/{{sample}}_genes.gff",
        faa=f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}/{{sample}}_proteins.faa",
        fna=f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}/{{sample}}_orfs.fna"

    conda:
        "../envs/prodigal.yaml"

    shell:"prodigal -i {input.contigs} -o {output.gff} -a {output.faa} -d {output.fna} -p meta -f gff"

rule whokaryote:
    input:
        contigs=f"{outdir}/results/03_assembly/size_filtered/{{sample}}_{minsize}/contigs_{{sample}}_{minsize}.fasta",
        prodigal_gff=f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}/{{sample}}_genes.gff"

    output:
        headers_euk=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/eukaryote_contig_headers.txt",
        headers_prok=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/prokaryote_contig_headers.txt",
        contigs_size=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/contigs_{minsize}.fasta",
        euk_fasta=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/eukaryotes.fasta",
        prok_fasta=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/prokaryotes.fasta",
        unclassified_fasta=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/unclassified.fasta"

    conda:
        "../envs/whokaryote.yaml"

    params:
        outdir=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}",
        minsize=config['min_contig_length']

    shell:
        "whokaryote.py --contigs {input.contigs} --outdir {params.outdir} --prodigal_file {input.prodigal_gff} --minsize {params.minsize} --f"


rule augustify:
    input:
        f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/eukaryotes.fasta"
    output:
        tax=f"{outdir}/results/04_gene_prediction/augustify/{{sample}}/{{sample}}_eukclass.txt",
        gff=f"{outdir}/results/04_gene_prediction/augustify/{{sample}}/{{sample}}_eukproteins.gff",
    conda:
        "../envs/augustus.yaml"
    params:
        param_file=os.path.abspath("resources/augustify_params.txt"),
        threads=config['threads'],
        script=os.path.abspath("workflow/scripts/augustify.py")
    shell:
        "python {params.script} -g {input} -p {params.param_file} -m {output.tax} -P {output.gff} -t {params.threads}"


