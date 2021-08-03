rule merge_per_treatment:
    input:
        forward=lambda wildcards: expand("scratch/host_filtering/{sample}_R1.fastq", project=config["project"], sample=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("scratch/filter/{sample}_R1.fasta", project=config["project"], sample=config["treatment"][wildcards.treatment]),
        rev=lambda wildcards: expand("scratch/host_filtering/{sample}_R2.fastq", project=config["project"], sample=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("scratch/filter/{sample}_R2.fasta", project=config["project"], sample=config["treatment"][wildcards.treatment]),
    output:
        forward = protected("scratch/treatment/{treatment}_forward.fastq"),
        rev = protected("scratch/treatment/{treatment}_rev.fastq"),
    log: 
        forward = "logs/pool/merge_per_treatment_{treatment}_forward.log",
        rev = "logs/pool/merge_per_treatment_{treatment}_rev.log"
    run: 
        shell("cat {input.forward}  > {output.forward} 2> {log.forward}")
        shell("cat {input.rev}  > {output.rev} 2> {log.rev}")
