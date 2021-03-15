rule predict_genes:
    input:
        "scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output:
        gbk="scratch/prodigal/coords.gbk",
        proteins="scratch/prodigal/proteins.faa",
        nucleotide="scratch/prodigal/orfs.fna",
    log:
        "logs/prodigal.log"
    conda: 
        "../../../envs/prodigal.yaml"
    shell: "prodigal -i {input} -o {output.gbk} -d {output.nucleotide} -a {output.proteins} -p meta"
