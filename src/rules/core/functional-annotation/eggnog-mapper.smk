rule move_proteins:
    input: "scratch/prodigal/proteins.faa"
    output: "results/annotation/emapper/proteins.faa"
    params:
        outdir="results/annotation/emapper"
    shell: "cp {input} {params.outdir}"

rule emapper_diamond:
    input:
        "results/annotation/emapper/proteins.faa"
    output:
        "results/annotation/emapper/proteins.faa.emapper.seed_orthologs"
    params:
        db=config["emapper_diamond"]
    conda:
        "../../../envs/eggnog-mapper.yaml"
    threads: 0
    shell: "emapper.py --dmnd_db {params.db} -m diamond --no_annot --no_file_comments --cpu {threads} -i {input} -o {input}"

rule emapper_annotation:
    input:
        seq="results/annotation/emapper/proteins.faa",
        diamond="results/annotation/emapper/proteins.faa.emapper.seed_orthologs"
    output:
        "results/annotation/emapper/proteins.faa.emapper.annotations"
    params:
        db=config['emapper_database']
    conda:
        "../../../envs/eggnog-mapper.yaml"
    threads: 0
    shell:
        """
        emapper.py --annotate_hits_table {input.diamond} --no_file_comments --cpu {threads} --data_dir {params.db} -o {input.seq}
        """
