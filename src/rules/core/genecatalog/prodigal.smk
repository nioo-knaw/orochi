"""
rule predict_genes:
    input:
        "scratch/assembly/{assembler}/minimus2/secondary.contigs.fasta"
    output:
        gbk="scratch/genecatalog/coords.gbk",
        proteins="scratch/genecatalog/proteins.faa",
        nucleotide="scratch/genecatalog/orfs.fna",
    log:
        "prodigal.log"
    conda: 
        "../../../envs/prodigal.yaml"
    shell: "prodigal -d {output.nucleotide} -a {output.proteins} -i {input} -m -o {output.gbk} -p anon"
"""
