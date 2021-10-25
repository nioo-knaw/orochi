rule predict_genes:
    input:
        "scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output:
        gbk=temp("scratch/prodigal/coords.gbk"),
        proteins=temp("scratch/prodigal/proteins.faa"),
        nucleotide=temp("scratch/prodigal/orfs.fna"),
    log:
        "logs/prodigal.log"
    conda: 
        "../../../envs/prodigal.yaml"
    shell: "prodigal -i {input} -o {output.gbk} -d {output.nucleotide} -a {output.proteins} -p meta"
