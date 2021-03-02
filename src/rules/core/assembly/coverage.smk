rule bamfiles:
    input:
        forward="scratch/host_filtering/{sample}_R1.fastq" if config['host_removal'] \
             else "scratch/filter/{sample}_R1.fasta",
        reverse="scratch/host_filtering/{sample}_R2.fastq" if config['host_removal'] \
             else "scratch/filter/{sample}_R2.fasta",
        assembly = "scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output:
        "scratch/coverm/bamfiles/secondary.contigs.fasta.{sample}_R1.fastq.bam"
    params:
        outdir="scratch/coverm/bamfiles"
    conda:
        "../../../envs/coverm.yaml"
    threads: 80
    shell:
        "coverm make -p bwa-mem -r {input.assembly} -1 {input.forward} -2 {input.reverse} -o {params.outdir} -t {threads}"

rule sort_readname:
    input:
        "scratch/coverm/bamfiles/secondary.contigs.fasta.{sample}_R1.fastq.bam"
    output:
        "scratch/coverm/bamfiles/readsorted/{sample}.bam"
    conda:
        "../../../envs/samtools.yaml"
    shell:
        "samtools sort -n {input} -o {output}"

rule coverage:
    input:
        expand("scratch/coverm/bamfiles/secondary.contigs.fasta.{sample}_R1.fastq.bam", sample=config["data"])
    output:
        protected("results/stats/coverage/coverage.tsv")
    conda:
        "../../../envs/coverm.yaml"
    threads: 80
    shell:
        "coverm contig -b {input} -t {threads} -o {output}"

#TO DO: Add stderr log?

"""
rule coverage_old:
    input:
        forward = expand("scratch/host_filtering/{sample}_R1.fastq", sample=config["data"]),
        reverse = expand("scratch/host_filtering/{sample}_R2.fastq", sample=config["data"]),
#        assembly = "scratch/assembly/megahit/all/meta-large/final.contigs.fa"
        # TODO: Decide what is the input here
        assembly=expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        table=protected("results/stats/coverage.tsv")
    conda:
        "../../../envs/coverm.yaml"
    threads: 80
    # TODO: Add log file stderr
    shell:
        "coverm contig --mapper bwa-mem --methods mean --reference {input.assembly} -1 {input.forward} -2 {input.reverse} --threads {threads} > {output.table}"
"""
