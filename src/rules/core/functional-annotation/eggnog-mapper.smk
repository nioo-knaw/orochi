rule eggnog_mapper_diamond:
    input:
        "scratch/genecatalog/{assembler}/{kmers}/allgenecalled.faa.gz"
    output:
        "scratch/genecatalog/{assembler}/{kmers}/allgenecalled.faa.gz.emapper.seed_orthologs"
    conda:
        "../../../envs/eggnog-mapper.yaml"
    threads: 16
    shell: "emapper.py --dmnd_db /data/db/eggnogdb/5.0.0/eggnog_proteins.dmnd -m diamond --no_annot --no_file_comments --cpu {threads} -i {input} -o {input}" 

rule eggnog_mapper_annotation:
    input:
        seq="scratch/genecatalog/{assembler}/{kmers}/allgenecalled.faa.gz",
        diamond="scratch/genecatalog/{assembler}/{kmers}/allgenecalled.faa.gz.emapper.seed_orthologs"
    output:
        "scratch/genecatalog/{assembler}/{kmers}/allgenecalled.faa.gz.emapper.annotations"
    conda:
        "../../../envs/eggnog-mapper.yaml"
    threads: 16
    shell: "emapper.py --annotate_hits_table {input.diamond} --no_file_comments --cpu {threads} --data_dir /scratch -o {input.seq}"