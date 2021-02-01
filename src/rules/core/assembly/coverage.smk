rule coverage:
    input:
        forward = expand("scratch/host_filtering/{sample}_R1.fastq", sample=config["data"]),
        reverse = expand("scratch/host_filtering/{sample}_R2.fastq", sample=config["data"]),
#        assembly = "scratch/assembly/megahit/all/meta-large/final.contigs.fa"
        # TODO: Decide what is the input here
        assembly=expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        "results/stats/coverage.tsv"
        #"scratch/coverm/coverage.tsv"
    params:
        bamdir="scratch/coverm/bamfiles"
    conda:
        "../../../envs/coverm.yaml"
    threads: 16
    shell:
        #"coverm contig --mapper bwa-mem --methods mean -c {input.forward} {input.reverse} --reference {input.assembly} -t {threads} --bam-file-cache-directory {params.bamdir}"
        "coverm contig -c {input.forward} {input.reverse} -r {input.assembly} -t {threads}"

#TO DO: Add stderr log?

rule bamfiles:
    input:
        forward = expand("scratch/host_filtering/{sample}_R1.fastq", sample=config["data"]),
        reverse = expand("scratch/host_filtering/{sample}_R2.fastq", sample=config["data"]),
        assembly = expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        "scratch/coverm/bamfiles/{sample}.assembly.fa.bam"
    params:
        outdir="scratch/coverm/bamfiles"
    conda:
        "../../../envs/coverm.yaml"
    threads: 16
    shell:
        "coverm contig -c {input.forward} {input.reverse} -r {input.assembly} --bam-file-cache-directory {params.outdir} -t {threads}"
