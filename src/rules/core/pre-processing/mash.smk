rule mash_sketch:
    input:
        forward="scratch/host_filtering/{sample}_R1.fastq" if config['host_removal'] \
            else "scratch/filter/{sample}_R1.fq",
        rev="scratch/host_filtering/{sample}_R2.fastq" if config['host_removal'] \
            else "scratch/filter/{sample}_R2.fq",
    output:
        temp("scratch/treatment/mash/{sample}.msh")
    params:
       prefix=lambda wildcards, output: os.path.join(os.path.dirname(str(output)),os.path.basename(str(output)))
    log: "logs/mash/mash_sketch_{sample}.log"
    conda: "../../../envs/mash.yaml"
    shell: "mash sketch -k 27 -s 10000 -o {params.prefix} -r {input.forward} {input.rev} 2> {log}"

rule mash_paste:
    input:
        expand("scratch/treatment/mash/{sample}.msh", sample=config["data"])
    output:
        temp("scratch/treatment/all.msh")
    params:
       prefix=lambda wildcards, output: os.path.join(os.path.dirname(str(output)), "all")
    log: "logs/mash/mash_paste.log"
    conda: "../../../envs/mash.yaml"
    shell: "mash paste {params.prefix} {input} 2> {log}"

rule mash_dist:
    input: "scratch/treatment/all.msh"
    output: temp("results/stats/mash/table.txt")
    log: "logs/mash/mash_dist.log"
    conda: "../../../envs/mash.yaml"
    shell: "mash dist -t {input} {input} > {output} 2> {log}"

rule draw_tree:
    input:
        txt="results/stats/mash/table.txt"
    output:
        svg="results/stats/mash/tree.svg"
    params:
        req="h_vmem=2G",
        samples=config["data"]
    log: "logs/mash/draw_tree.log"
    conda:
        "../../../envs/tree.yaml"
    script:
        "../../../scripts/draw_tree.py"

