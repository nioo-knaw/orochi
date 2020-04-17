rule merge_per_treatment:
    input:
#        forward=lambda wildcards: expand("scratch/filter/{sample}_forward_paired.fq.gz", project=config["project"], sample=config["treatment"][wildcards.treatment])
        forward=lambda wildcards: expand("scratch/host_filtering/{sample}_R1.fastq", project=config["project"], sample=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("scratch/filter/{sample}_R1.fasta", project=config["project"], sample=config["treatment"][wildcards.treatment]),
        reverse=lambda wildcards: expand("scratch/host_filtering/{sample}_R2.fastq", project=config["project"], sample=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("scratch/filter/{sample}_R2.fasta", project=config["project"], sample=config["treatment"][wildcards.treatment]),
#        unpaired=lambda wildcards: expand("scratch/host_filtering/{sample}_unpaired_filtered.fastq", project=config["project"], sample=config["treatment"][wildcards.treatment]) if config['host_removal'] \
#             else expand("scratch/filter/{sample}_forward_unpaired.fq.gz", project=config["project"], sample=config["treatment"][wildcards.treatment])
    output:
        forward = protected("scratch/treatment/{treatment}_forward.fastq"),
        reverse = protected("scratch/treatment/{treatment}_reverse.fastq"),
#        unpaired = protected("scratch/treatment/{treatment}_unpaired.fastq")
    run: 
        if config['host_removal']:
            shell("cat {input.forward}  > {output.forward}")
            shell("cat {input.reverse}  > {output.reverse}")
#            shell("cat {input.unpaired} > {output.unpaired}")
        else:
            shell("cat {input.forward}  > {output.forward}")
            shell("cat {input.reverse}  > {output.reverse}")
#            shell("cat {input.unpaired} > {output.unpaired}")

