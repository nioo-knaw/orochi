rule merge_per_treatment:
    input:
        forward=lambda wildcards: expand("scratch/host_filtering/{sample}_R1.fastq", project=config["project"], sample=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("scratch/filter/{sample}_R1.fasta", project=config["project"], sample=config["treatment"][wildcards.treatment]),
        reverse=lambda wildcards: expand("scratch/host_filtering/{sample}_R2.fastq", project=config["project"], sample=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("scratch/filter/{sample}_R2.fasta", project=config["project"], sample=config["treatment"][wildcards.treatment]),
    output:
        forward = protected("scratch/treatment/{treatment}_forward.fastq"),
        reverse = protected("scratch/treatment/{treatment}_reverse.fastq"),
    run: 
        if config['host_removal']:
            shell("cat {input.forward}  > {output.forward}")
            shell("cat {input.reverse}  > {output.reverse}")
        else:
            shell("cat {input.forward}  > {output.forward}")
            shell("cat {input.reverse}  > {output.reverse}")

rule sort_per_treatment:
    input:
        forward=lambda wildcards: expand("scratch/host_filtering/{sample}_R1.fastq", sample=config["data"], treatment=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("scratch/filter/{sample}_R1.fasta", sample=config["data"], treatment=config["treatment"][wildcards.treatment]),
        reverse=lambda wildcards: expand("scratch/host_filtering/{sample}_R2.fastq", sample=config["data"], treatment=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("scratch/filter/{sample}_R2.fasta", sample=config["data"], treatment=config["treatment"][wildcards.treatment]),
    output:
        forward = protected("scratch/sort/{treatment}/{sample}_R1.fastq"),
        reverse = protected("scratch/sort/{treatment}/{sample}_R2.fastq"),
    params:
        outdir="scratch/sort/{treatment}"
    run: 
        if config['host_removal']:
            shell("cp {input.forward} {params.outdir}")
            shell("cp {input.reverse} {params.outdir}")
        else:
            shell("cp {input.forward} {params.outdir}")
            shell("cp {input.reverse} {params.outdir}")
