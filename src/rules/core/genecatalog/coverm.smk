rule coverage:
    input:
        forward = expand("scratch/host_filtering/{sample}_R1.fastq", sample=config["data"]),
        reverse = expand("scratch/host_filtering/{sample}_R2.fastq", sample=config["data"]),
#        assembly = "scratch/assembly/megahit/all/meta-large/final.contigs.fa"
        # TODO: Decide what is the input here
        assembly=expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        table="results/stats/coverage.tsv"
        bam="scratch/coverm/{assembler}/{treatment}/{kmers}/assembly.fa.{sample}_1.fastq.bam"
    conda:
        "../../../envs/coverm.yaml"
    threads: 16
    # TODO: Add log file stderr
    shell:
        "coverm contig --mapper bwa-mem --methods mean --bam-file-cache-directory --reference {input.assembly} -1 {input.forward} -2 {input.reverse} --threads {threads} > {output.table}"

"""
rule coverm_treatment:
    input:
        contigs="scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa", 
        forward = "scratch/treatment/{treatment}_forward.fastq",
        reverse = "scratch/treatment/{treatment}_reverse.fastq",
        index="scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.bwt"
    output:
         "scratch/coverm/{assembler}/{treatment}/{kmers}/assembly.fa.{treatment}.fastq.bam",
    log:
        "scratch/coverm/{assembler}/{treatment}/{kmers}/{treatment}.log"
    params:
        outdir="scratch/coverm/{assembler}/{treatment}/{kmers}"
    threads: 16
    conda:
        "../../../envs/coverm.yaml"
    shell: "coverm make -r {input.contigs} -c {input.forward} {input.reverse} -o {params.outdir} -t {threads} 2> {log}"
"""
"""
rule relocate_sample:
    input:
        forward = expand("scratch/unpack/{sample}_1.fastq", sample=config["data"]),
        reverse = expand("scratch/unpack/{sample}_2.fastq", sample=config["data"])
    output:
        forward = "scratch/assembly/{assembler}/{treatment}/{kmers}/{sample}_1.fastq",
        reverse = "scratch/assembly/{assembler}/{treatment}/{kmers}/{sample}_2.fastq"
    threads: 16
    run:
        shell("cp {input.forward} {output.forward}")
        shell("cp {input.reverse} {output.reverse}")              
"""
"""
rule coverm_sample:
    input:
        contigs="scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa", 
        forward = expand("scratch/unpack/{sample}_1.fastq", sample=config["data"]),
        reverse = expand("scratch/unpack/{sample}_2.fastq", sample=config["data"]),
        index="scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.bwt"
    output:
        "scratch/coverm/{assembler}/{treatment}/{kmers}/coveragebysample.txt"
    log:
        "scratch/coverm/{assembler}/{treatment}/{kmers}/coverm_sample.log"
    params:
        outdir="scratch/coverm/{assembler}/{treatment}/{kmers}"
    threads: 16
    conda:
        "../../../envs/coverm.yaml"
    shell: "coverm contig --methods count --mapper minimap2-sr --proper-pairs-only -1 {input.forward} -2 {input.reverse} --reference {input.contigs} --threads {threads} 2> {log} > {output}"
"""
