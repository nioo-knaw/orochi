# Trim adapters and low quality regions
rule trimmomatic:
    input:
        r1="{project}/unpack/{sample}_1.fastq.gz",
        r2="{project}/unpack/{sample}_2.fastq.gz",
    output:
        r1=protected("{project}/trimmomatic/{sample}_forward_paired.fq.gz"),
        r1_unpaired=protected("{project}/trimmomatic/{sample}_forward_unpaired.fq.gz"),
        r2=protected("{project}/trimmomatic/{sample}_reverse_paired.fq.gz"),
        r2_unpaired=protected("{project}/trimmomatic/{sample}_reverse_unpaired.fq.gz"),
    log:
        "{project}/trimmomatic/{sample}.log"
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
        fw_unpaired=temp("{project}/trimmomatic/{sample}_forward_unpaired.fq.gz"),
        rev_unpaired=temp("{project}/trimmomatic/{sample}_reverse_unpaired.fq.gz"),
    output:
        "{project}/trimmomatic/{sample}_unpaired_combined.fq.gz"
    shell: "zcat {input} | gzip -c > {output}"

rule merge_per_treatment:
    input:
#        forward=lambda wildcards: expand("{project}/trimmomatic/{sample}_forward_paired.fq.gz", project=config["project"], sample=config["treatment"][wildcards.treatment])
        forward=lambda wildcards: expand("{project}/host_filtering/{sample}_R1_paired_filtered.fastq", project=config["project"], sample=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("{project}/trimmomatic/{sample}_forward_paired.fq.gz", project=config["project"], sample=config["treatment"][wildcards.treatment]),
        reverse=lambda wildcards: expand("{project}/host_filtering/{sample}_R2_paired_filtered.fastq", project=config["project"], sample=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("{project}/trimmomatic/{sample}_reverse_paired.fq.gz", project=config["project"], sample=config["treatment"][wildcards.treatment]),
        unpaired=lambda wildcards: expand("{project}/host_filtering/{sample}_unpaired_filtered.fastq", project=config["project"], sample=config["treatment"][wildcards.treatment]) if config['host_removal'] \
             else expand("{project}/trimmomatic/{sample}_forward_unpaired.fq.gz", project=config["project"], sample=config["treatment"][wildcards.treatment])
    output:
        forward = temporary("{project}/treatment/{treatment}_forward.fastq"),
        reverse = temporary("{project}/treatment/{treatment}_reverse.fastq"),
        unpaired = temporary("{project}/treatment/{treatment}_unpaired.fastq")
    run: 
        if config['host_removal']:
            shell("cat {input.forward}  > {output.forward}")
            shell("cat {input.reverse}  > {output.reverse}")
            shell("cat {input.unpaired} > {output.unpaired}")
        else:
            shell("cat {input.forward}  > {output.forward}")
            shell("cat {input.reverse}  > {output.reverse}")
            shell("cat {input.unpaired} > {output.unpaired}")

