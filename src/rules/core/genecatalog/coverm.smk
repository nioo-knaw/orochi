rule coverage:
    input:
        forward = expand("scratch/host_filtering/{sample}_R1.fastq", sample=config["data"]),
        reverse = expand("scratch/host_filtering/{sample}_R2.fastq", sample=config["data"]),
#        assembly = "scratch/assembly/megahit/all/meta-large/final.contigs.fa"
        # TODO: Decide what is the input here
        assembly=expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        "scratch/genecatalog/coverage.tsv"
    conda:
        "../../../envs/coverm.yaml"
    threads: 16
    # TODO: Add log file stderr
    shell:
        "coverm contig --mapper bwa-mem --methods mean --reference {input.assembly} -1 {input.forward} -2 {input.reverse} --threads {threads} > {output}"

rule coverm_treatment:
    input:
        contigs="scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa", 
        forward = "scratch/treatment/{treatment}_forward.fastq",
        reverse = "scratch/treatment/{treatment}_reverse.fastq",
        index="scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.bwt"
    output:
         "scratch/coverm/{assembler}/{treatment}/{kmers}/assembly.{treatment}_forwardstq.bam",
         "scratch/coverm/{assembler}/{treatment}/{kmers}/assembly.{treatment}_forwardstq.bam.bai"
    log:
        "scratch/coverm/{assembler}/{treatment}/{kmers}/{treatment}.log"
    params:
        outdir="scratch/coverm/{assembler}/{treatment}/{kmers}"
    threads: 16
    conda:
        "../../../envs/coverm.yaml"
    shell: "coverm make -r {input.contigs} -c {input.forward} {input.reverse} -o {params.outdir} -t {threads} 2> {log}"