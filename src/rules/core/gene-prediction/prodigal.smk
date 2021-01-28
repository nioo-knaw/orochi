rule predict_genes:
    input:
        #"scratch/assembly/megahit/minimus2/secondary.contigs.fasta" if config["multiple-assemblies"] else 
        expand("scratch/assembly/megahit/{treatment}/{kmers}/final.contigs.fa", treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        gbk="scratch/genecatalog/coords.gbk",
        proteins="scratch/genecatalog/proteins.faa",
        nucleotide="scratch/genecatalog/orfs.fna",
    log:
        "prodigal.log"
    conda: 
        "../../../envs/prodigal.yaml"
    shell: "prodigal -i {input} -o {output.gbk} -d {output.nucleotide} -a {output.proteins} -p meta"
