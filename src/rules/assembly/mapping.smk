rule treatment_bwa_index:
    input:
         "scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa"
    output:
        "scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.bwt"
    log: "scratch/assembly/{assembler}/{treatment}/{kmers}/bwa-index.log"
    shell: "/data/tools/bwa/default/bin/bwa index {input} > {log}"

rule bamm_treatment:
    input:
        contigs="scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa", 
        forward = "scratch/treatment/{treatment}_forward.fastq",
        reverse = "scratch/treatment/{treatment}_reverse.fastq",
        unpaired = "scratch/treatment/{treatment}_unpaired.fastq",
        index="scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.bwt"
    output:
         "scratch/bamm/{assembler}/{treatment}/{kmers}/assembly.{treatment}_forwardstq.bam",
         "scratch/bamm/{assembler}/{treatment}/{kmers}/assembly.{treatment}_forwardstq.bam.bai"
    log:
        "scratch/bamm/{assembler}/{treatment}/{kmers}/{treatment}.log"
    params:
        outdir="scratch/bamm/{assembler}/{treatment}/{kmers}"
    threads: 16
    conda:
        "../../envs/bamm.yaml"
    shell: "bamm make --keep_unmapped --kept -d {input.contigs} -c {input.forward} {input.reverse} -s {input.unpaired} -o {params.outdir} -t {threads} 2> {log}"


rule bamm_samples:
    input:
        contigs="scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa",
        index="scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.bwt",
        forward = "scratch/host_filtering/{sample}_R1_paired_filtered.fastq" if config['host_removal'] else \
        "scratch/trimmomatic/{sample}_forward_paired.fq.gz",
        reverse = "scratch/host_filtering/{sample}_R2_paired_filtered.fastq" if config['host_removal'] else \
       "scratch/trimmomatic/{sample}_reverse_paired.fq.gz",
        unpaired = "scratch/host_filtering/{sample}_unpaired_filtered.fastq" if config['host_removal'] else "scratch/trimmomatic/{sample}_unpaired_combined.fq.gz",
    output:
        "scratch/bamm/{assembler}/{treatment}/{kmers}/assembly.{sample}_R1_paired_filteredstq.bam" if config['host_removal'] else "scratch/bamm/{assembler}/{treatment}/{kmers}/assembly.{sample}_forward_paired.bam", 
        "scratch/bamm/{assembler}/{treatment}/{kmers}/assembly.{sample}_R1_paired_filteredstq.bam.bai" if config['host_removal'] else "scratch/bamm/{assembler}/{treatment}/{kmers}/assembly.{sample}_forward_paired.bam.bai", 
    log:
        "scratch/bamm/{sample}_{assembler}_{treatment}_{kmers}.log"
    params:
        outdir="scratch/bamm/{assembler}/{treatment}/{kmers}/"
    threads: 16
    conda:
        "../../envs/bamm.yaml"
    shell: "bamm make --keep_unmapped --kept -d {input.contigs} -c {input.forward} {input.reverse} -s {input.unpaired} -o {params.outdir} -t {threads} 2> {log}"

