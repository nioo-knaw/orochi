rule coverage:
    input:
        forward=lambda wildcards: expand("scratch/host_filtering/{sample}_R1.fastq", sample=config["treatment"][wildcards.treatment]),
        reverse=lambda wildcards: expand("scratch/host_filtering/{sample}_R2.fastq", sample=config["treatment"][wildcards.treatment]),
#        assembly = "scratch/assembly/megahit/all/meta-large/final.contigs.fa"
        # TODO: Decide what is the input here
        assembly=expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        "results/stats/coverage.tsv"
        #"scratch/coverm/coverage.tsv"
    conda:
        "../../../envs/coverm.yaml"
    threads: 16
    shell:
        #"coverm contig --mapper bwa-mem --methods mean -c {input.forward} {input.reverse} --reference {input.assembly} -t {threads} --bam-file-cache-directory {params.bamdir}"
        "coverm contig --mapper bwa-mem -c {input.forward} {input.reverse} -r {input.assembly} -t {threads} -o {output}"

#TO DO: Add stderr log?

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
        "coverm contig -c {input.forward} {input.reverse} -r {input.assembly} --bam-file-cache-directory {params.outdir} -t {threads}"
