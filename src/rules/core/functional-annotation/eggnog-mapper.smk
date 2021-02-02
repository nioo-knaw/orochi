"""
rule fetch_emapper_data:
    input:  "download_eggnog_data.py"
    output: "data/db/eggnogdb/5.0.0/eggnog_proteins.dmnd"
    conda:
        "../../../envs/eggnog-mapper.yaml"
    shell: "python download_eggnog_data.py"
"""

rule eggnog_mapper_diamond:
    input:
        "scratch/genecatalog/proteins.faa"
    output:
        "scratch/genecatalog/proteins.faa.emapper.seed_orthologs"
    conda:
        "../../../envs/eggnog-mapper.yaml"
    threads: 16
    shell: 
        "python download_eggnog_data.py | emapper.py --dmnd_db /data/db/eggnogdb/5.0.0/eggnog_proteins.dmnd -m diamond --no_annot --no_file_comments --cpu {threads} -i {input} -o {input}"

rule eggnog_mapper_annotation:
    input:
        seq="scratch/genecatalog/proteins.faa",
        diamond="scratch/genecatalog/proteins.faa.emapper.seed_orthologs"
    output:
        "results/annotation/eggnog-mapper/proteins.faa.emapper.annotations"
    conda:
        "../../../envs/eggnog-mapper.yaml"
    threads: 16
    shell: "emapper.py --annotate_hits_table {input.diamond} --no_file_comments --cpu {threads} --data_dir /scratch -o {input.seq}"
