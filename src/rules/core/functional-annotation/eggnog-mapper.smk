rule eggnog_mapper_diamond:
    input:
        "scratch/genecatalog/proteins.faa"
    output:
        "scratch/genecatalog/proteins.faa.emapper.seed_orthologs"
    conda:
        "../../../envs/eggnog-mapper.yaml"
    threads: 16
    shell: "emapper.py --dmnd_db /data/db/eggnogdb/5.0.0/eggnog_proteins.dmnd -m diamond --no_annot --no_file_comments --cpu {threads} -i {input} -o {input}"

rule eggnog_mapper_annotation:
    input:
        seq="scratch/genecatalog/proteins.faa",
        diamond="scratch/genecatalog/proteins.faa.emapper.seed_orthologs"
    output:
        "results/annotation/eggnog-mapper/proteins.faa.emapper.annotations"
    params:
        datadir="data/db/eggnogdb/5.0.0"
    conda:
        "../../../envs/eggnog-mapper.yaml"
    threads: 16
    shell:
        """
        download_eggnog_data.py --data_dir {params.datadir} -y
        emapper.py --annotate_hits_table {input.diamond} --no_file_comments --cpu {threads} --data_dir {params.datadir} -o {input.seq}
        """
#Having the input and output directory differ here may cause problems because of how emapper names output using input files as prefix.
