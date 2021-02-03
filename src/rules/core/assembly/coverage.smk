rule bamfiles:
    input:
        forward=lambda wildcards: expand("scratch/host_filtering/{sample}_R1.fastq", sample=config["treatment"][wildcards.treatment]),
        reverse=lambda wildcards: expand("scratch/host_filtering/{sample}_R2.fastq", sample=config["treatment"][wildcards.treatment]),
        assembly = expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa", treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        forward = "scratch/coverm/bamfiles/{treatment}/assembly.fa.{sample}_R1.fastq.bam",
        reverse = "scratch/coverm/bamfiles/{treatment}/assembly.fa.{sample}_R2.fastq.bam"
    params:
        outdir="scratch/coverm/bamfiles/{treatment}"
    conda:
        "../../../envs/coverm.yaml"
    threads: 16
    shell:
        "coverm make -p bwa_mem -r {input.assembly} -c {input.forward} {input.reverse} -o {params.outdir} -t {threads}"

rule coverage:
    input:
        forward = expand("scratch/coverm/bamfiles/{treatment}/assembly.fa.{sample}_R1.fastq.bam", treatment=config["treatment"], sample=config["data"]),
        reverse = expand("scratch/coverm/bamfiles/{treatment}/assembly.fa.{sample}_R2.fastq.bam", treatment=config["treatment"], sample=config["data"]),
        assembly=expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        "results/stats/coverage.tsv"
    conda:
        "../../../envs/coverm.yaml"
    threads: 16
    shell:
        "coverm contig -b {input.forward} {input.reverse} -r {input.assembly} -t {threads} -o {output}"

#TO DO: Add stderr log?
