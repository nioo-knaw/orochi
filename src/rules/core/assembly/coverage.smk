rule bamfiles:
    input:
        forward = "scratch/host_filtering/{sample}_R1.fastq",
        reverse = "scratch/host_filtering/{sample}_R2.fastq",
        assembly = "scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output:
        "scratch/coverm/bamfiles/secondary.contigs.fasta.{sample}_R1.fastq.bam"
    params:
        outdir="scratch/coverm/bamfiles/"
    conda:
        "../../../envs/coverm.yaml"
    threads: 16
    shell:
        "coverm make -p bwa-mem -r {input.assembly} -1 {input.forward} -2 {input.reverse} -o {params.outdir} -t {threads}"

rule coverage:
    input: expand("scratch/coverm/bamfiles/secondary.contigs.fasta.{sample}_R1.fastq.bam", sample=config["data"])
    output:
        "results/stats/coverage.tsv"
    conda:
        "../../../envs/coverm.yaml"
    threads: 16
    shell:
        "coverm contig -b {input} -t {threads} -o {output}"

#TO DO: Add stderr log?
