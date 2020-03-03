# Trim adapters and low quality regions
rule trimmomatic:
    input:
        r1="scratch/unpack/{sample}_1.fastq.gz",
        r2="scratch/unpack/{sample}_2.fastq.gz",
    output:
        r1=protected("scratch/trimmomatic/{sample}_forward_paired.fq.gz"),
        r1_unpaired=protected("scratch/trimmomatic/{sample}_forward_unpaired.fq.gz"),
        r2=protected("scratch/trimmomatic/{sample}_reverse_paired.fq.gz"),
        r2_unpaired=protected("scratch/trimmomatic/{sample}_reverse_unpaired.fq.gz"),
    log:
        "scratch/trimmomatic/{sample}.log"
    threads: 16 
    params:
        # list of trimmers (see manual)
        #trimmer=["ILLUMINACLIP:{config["adapters"]}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:100"],
        # optional parameters
        extra="-threads {snakemake.threads}"
    wrapper:
        "0.22.0/bio/trimmomatic/pe"


rule trimmomatic_combine_unpaired:
    input:
        fw_unpaired="scratch/trimmomatic/{sample}_forward_unpaired.fq.gz",
        rev_unpaired="scratch/trimmomatic/{sample}_reverse_unpaired.fq.gz",
    output:
        "scratch/trimmomatic/{sample}_unpaired_combined.fq.gz"
    shell: "zcat {input} | gzip -c > {output}"

rule merge_per_treatment:
    input:
#        forward=lambda wildcards: expand("scratch/trimmomatic/{sample}_forward_paired.fq.gz", project=config["project"], sample=config["treatment"][wildcards.treatment])
        forward=lambda wildcards: expand("scratch/host_filtering/{sample}_R1_paired_filtered.fastq", project=config["project"], sample=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("scratch/trimmomatic/{sample}_forward_paired.fq.gz", project=config["project"], sample=config["treatment"][wildcards.treatment]),
        reverse=lambda wildcards: expand("scratch/host_filtering/{sample}_R2_paired_filtered.fastq", project=config["project"], sample=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("scratch/trimmomatic/{sample}_reverse_paired.fq.gz", project=config["project"], sample=config["treatment"][wildcards.treatment]),
        unpaired=lambda wildcards: expand("scratch/host_filtering/{sample}_unpaired_filtered.fastq", project=config["project"], sample=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("scratch/trimmomatic/{sample}_forward_unpaired.fq.gz", project=config["project"], sample=config["treatment"][wildcards.treatment])
    output:
        forward = protected("scratch/treatment/{treatment}_forward.fastq"),
        reverse = protected("scratch/treatment/{treatment}_reverse.fastq"),
        unpaired = protected("scratch/treatment/{treatment}_unpaired.fastq")
    run: 
        if config['host_removal']:
            shell("cat {input.forward}  > {output.forward}")
            shell("cat {input.reverse}  > {output.reverse}")
            shell("cat {input.unpaired} > {output.unpaired}")
        else:
            shell("cat {input.forward}  > {output.forward}")
            shell("cat {input.reverse}  > {output.reverse}")
            shell("cat {input.unpaired} > {output.unpaired}")

