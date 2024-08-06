""" Rules related to gene prediction """

rule prodigal:
    input:
        contigs=branch(config['assembly_method'] == "coassembly",
            then="results/03_assembly/coassembly/assembly_{sample}/{sample}_assembly.fasta",
            otherwise="results/03_assembly/single_sample_assembly/{sample}/{sample}_assembly.fasta")

    output:
        gff="results/04_gene_prediction/prodigal/{sample}/{sample}_genes.gff",
        faa="results/04_gene_prediction/prodigal/{sample}/{sample}_proteins.faa",
        fna="results/04_gene_prediction/prodigal/{sample}/{sample}_orfs.fna"

    conda:
        "../envs/prodigal.yaml"

    shell:"prodigal -i {input.contigs} -o {output.gff} -a {output.faa} -d {output.fna} -p meta -f gff"


rule downstream_test:
    input:
        "results/04_gene_prediction/prodigal/{sample}/{sample}_genes.gff"

    output:
        test_file="results/05_test/{sample}/{sample}_test.txt"
    run:
        shell("touch {output.test_file}")