rule treatment_bwa_index:
    input:
         "scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa"
    output:
        "scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.bwt"
    log: "scratch/assembly/{assembler}/{treatment}/{kmers}/bwa-index.log"
    shell: "/data/tools/bwa/default/bin/bwa index {input} > {log}"

"""
rule bamm_treatment:
    input:
        contigs="scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa", 
        forward = "scratch/treatment/{treatment}_forward.fastq",
        reverse = "scratch/treatment/{treatment}_reverse.fastq",
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
        "../../../envs/bamm.yaml"
    shell: "bamm make --keep_unmapped --kept -d {input.contigs} -c {input.forward} {input.reverse} -o {params.outdir} -t {threads} 2> {log}"

rule bamm_samples:
    input:
        contigs="scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa",
        index="scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.bwt",
        forward = "scratch/host_filtering/{sample}_R1.fastq" if config['host_removal'] else \
        "scratch/filter/{sample}_R1.fasta",
        reverse = "scratch/host_filtering/{sample}_R2.fastq" if config['host_removal'] else \
       "scratch/filter/{sample}_R2.fasta",
    output:
        "scratch/bamm/{assembler}/{treatment}/{kmers}/assembly.{sample}_R1_paired_filteredstq.bam" if config['host_removal'] else "scratch/bamm/{assembler}/{treatment}/{kmers}/assembly.{sample}_forward_paired.bam", 
        "scratch/bamm/{assembler}/{treatment}/{kmers}/assembly.{sample}_R1_paired_filteredstq.bam.bai" if config['host_removal'] else "scratch/bamm/{assembler}/{treatment}/{kmers}/assembly.{sample}_forward_paired.bam.bai", 
    log:
        "scratch/bamm/{sample}_{assembler}_{treatment}_{kmers}.log"
    params:
        outdir="scratch/bamm/{assembler}/{treatment}/{kmers}/"
    threads: 16
    conda:
        "../../../envs/bamm.yaml"
    shell: "bamm make --keep_unmapped --kept -d {input.contigs} -c {input.forward} {input.reverse} -o {params.outdir} -t {threads} 2> {log}"
"""