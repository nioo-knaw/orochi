rule rename_megahit:
    input:
        "{project}/assembly/megahit/{treatment}/{kmers}/final.contigs.fa"
    output:
        gzip="{project}/assembly/megahit/{treatment}/{kmers}/assembly.fa.gz",
        fasta=temp("{project}/assembly/megahit/{treatment}/{kmers}/assembly.fa")
    run:
        shell("zcat {input} | awk '{{print $1}}' | sed 's/_/contig/' > {output.fasta}")
        shell("gzip -c {output.fasta} > {output.gzip}")

rule rename_spades:
    input:
        "{project}/assembly/spades/{treatment}/{kmers}/contigs.fasta"
    output:
        gzip=protected("{project}/assembly/spades/{treatment}/{kmers}/assembly.fa.gz"),
        fasta=temp("{project}/assembly/spades/{treatment}/{kmers}/assembly.fa")
    run:
        shell("cat {input} | awk '{{print $1}}' | sed 's/_/contig/' > {output.fasta}")
        shell("gzip -c {output.fasta} > {output.gzip}")

rule rename_idba:
    input:
        "{project}/assembly/idba/scaffold.fa"
    output:
        gzip="{project}/assembly/idba/assembly.fa.gz",
        fasta=temp("{project}/assembly/idba/assembly.fa")
    run:
        shell("cat {input} | awk '{{print $1}}' | sed 's/_/contig/' > {output.fasta}")
        shell("gzip -c {output.fasta} > {output.gzip}")

rule treatment_bwa_index:
    input:
         "{project}/assembly/{assembler}/{treatment}/{kmers}/assembly.fa"
    output:
        "{project}/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.bwt"
    log: "{project}/assembly/{assembler}/{treatment}/{kmers}/bwa-index.log"
    shell: "/data/tools/bwa/default/bin/bwa index {input} > {log}"

rule bamm_treatment:
    input:
        contigs="{project}/assembly/{assembler}/{treatment}/{kmers}/assembly.fa", 
        forward = "{project}/treatment/{treatment}_forward.fastq",
        reverse = "{project}/treatment/{treatment}_reverse.fastq",
        unpaired = "{project}/treatment/{treatment}_unpaired.fastq",
        index="{project}/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.bwt"
    output:
         "{project}/bamm/{assembler}/{treatment}/{kmers}/assembly.{treatment}_forwardstq.bam",
         "{project}/bamm/{assembler}/{treatment}/{kmers}/assembly.{treatment}_forwardstq.bam.bai"
    log:
        "{project}/bamm/{assembler}/{treatment}/{kmers}/{treatment}.log"
    params:
        outdir="{project}/bamm/{assembler}/{treatment}/{kmers}"
    threads: 16
    conda:
        "../../envs/bamm.yaml"
    shell: "bamm make --keep_unmapped --kept -d {input.contigs} -c {input.forward} {input.reverse} -s {input.unpaired} -o {params.outdir} -t {threads} 2> {log}"


rule bamm_samples:
    input:
        contigs="{project}/assembly/{assembler}/{treatment}/{kmers}/assembly.fa",
        index="{project}/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.bwt",
        forward = "{project}/host_filtering/{sample}_R1_paired_filtered.fastq" if config['host_removal'] else \
        "{project}/trimmomatic/{sample}_forward_paired.fq.gz",
        reverse = "{project}/host_filtering/{sample}_R2_paired_filtered.fastq" if config['host_removal'] else \
       "{project}/trimmomatic/{sample}_reverse_paired.fq.gz",
        unpaired = "{project}/host_filtering/{sample}_unpaired_filtered.fastq" if config['host_removal'] else "{project}/trimmomatic/{sample}_unpaired_combined.fq.gz",
    output:
        "{project}/bamm/{assembler}/{treatment}/{kmers}/assembly.{sample}_R1_paired_filteredstq.bam" if config['host_removal'] else "{project}/bamm/{assembler}/{treatment}/{kmers}/assembly.{sample}_forward_paired.bam", 
        "{project}/bamm/{assembler}/{treatment}/{kmers}/assembly.{sample}_R1_paired_filteredstq.bam.bai" if config['host_removal'] else "{project}/bamm/{assembler}/{treatment}/{kmers}/assembly.{sample}_forward_paired.bam.bai", 
    log:
        "{project}/bamm/{sample}.log"
    params:
        outdir="{project}/bamm/"
    threads: 16
    shell: "source /data/tools/samtools/1.3/env.sh; source /data/tools/BamM/1.7.0/env.sh; bamm make --keep_unmapped --kept -d {input.contigs} -c {input.forward} {input.reverse} -s {input.unpaired} -o {params.outdir} -t {threads} 2> {log}"

