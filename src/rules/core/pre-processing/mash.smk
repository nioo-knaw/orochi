rule mash_sketch:
    input:
        forward="scratch/host_filtering/{sample}_R1.fastq" if config['host_removal'] \
            else "scratch/filter/{sample}_R1.fasta",
        reverse="scratch/host_filtering/{sample}_R2.fastq" if config['host_removal'] \
            else "scratch/filter/{sample}_R2.fasta",
    output:
        "scratch/treatment/mash/{sample}.msh"
    params:
       prefix="scratch/treatment/mash/{sample}"
    conda: "../../../envs/mash.yaml"
    shell: "mash sketch -k 27 -s 10000 -o {params.prefix} -r {input.forward} {input.reverse}"

rule mash_paste:
    input:
        expand("scratch/treatment/mash/{sample}.msh", sample=config["data"])
    output:
        "scratch/treatment/all.msh"
    params:
       prefix="scratch/treatment/all"
    conda: "../../../envs/mash.yaml"
    shell: "mash paste {params.prefix} {input}"

rule mash_dist:
    input: "scratch/treatment/all.msh"
    output: "results/stats/mash/table.txt"
    conda: "../../../envs/mash.yaml"
    shell: "mash dist -t {input} {input} > {output}"

rule draw_tree:
    input:
        txt="results/stats/mash/table.txt"
    output:
        svg="results/stats/mash/tree.svg"
    params:
        req="h_vmem=2G",
        samples=config["data"]
    conda:
        "../../../envs/tree.yaml"
    script:
        "../../../scripts/draw_tree.py"

