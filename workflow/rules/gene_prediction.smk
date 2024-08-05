""" Rules related to gene prediction """

rule prodigal:
    input:
        contigs=branch(config['assembly_method'] == "coassembly",
            then="results/03_assembly/coassembly/assembly_{sample_pool}/{sample_pool}_assembly.fasta",
            otherwise="results/03_assembly/single_sample_assembly/{sample}/assembly.fasta"),


    output:
        gff="results/04_gene_prediction/prodigal/{sample_pool}/{sample_pool}_genes.gff",
        faa="results/04_gene_prediction/prodigal/{sample_pool}/{sample_pool}_proteins.faa",
        fna="results/04_gene_prediction/prodigal/{sample_pool}/{sample_pool}_orfs.fna"


    conda:
        "../envs/prodigal.yaml"

    shell:
        "prodigal -i input.contigs -o output.gff -a output.faa -d output.fna -p meta -f gff"



