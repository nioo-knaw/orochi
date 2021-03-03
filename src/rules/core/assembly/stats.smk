rule quast:
    input:
        "scratch/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.gz"
    output:
        quast="scratch/assembly/{assembler}/{treatment}/{kmers}/quast/report.txt",
    params:
        outdir="scratch/assembly/{assembler}/{treatment}/{kmers}/quast"
    log:
        "logs/assembly/{assembler}/{treatment}/{kmers}/quast/quast.log"
    conda:
        "../../../envs/quast.yaml"
    threads: 80
    shell: "metaquast.py -o {params.outdir} --min-contig 0 --max-ref-number 0 -t {threads} {input} 2>&1 > {log}"

rule quast_format:
    input:
        "scratch/assembly/{assembler}/{treatment}/{kmers}/quast/report.txt"
    output:
        "scratch/stats/{assembler}/{treatment}/{kmers}/quast.report.txt"
    log:
        "logs/stats/{assembler}/{treatment}/{kmers}/quast_format.log"
    params:
       run="{assembler}-{treatment}-{kmers}"
    shell: "printf '{params.run}\t' > {output} && cat {input} | sed 's/   */:/g' | cut -d : -f 2 | tr '\n' '\t' | cut -f 2- >> {output} 2> {log}"

rule quast_merge:
    input:
        quast = expand("scratch/stats/{assembler}/{treatment}/{kmers}/quast.report.txt", assembler=config["assembler"], treatment=config["treatment"], kmers=config["assembly-klist"]),
        full = expand("scratch/assembly/{assembler}/{treatment}/{kmers}/quast/report.txt", assembler=config["assembler"], treatment=config["treatment"], kmers=config["assembly-klist"]),
    output:
        protected("results/stats/quast.report.txt")
    log:
        "logs/stats/quast_merge.log"
    run:
         # Get only the first origin quast output file and get the first column for use as header
         firstfile = input.full[0]
         shell("cat {firstfile} | sed 's/   */:/g' | cut -d : -f 1 | tr '\n' '\t' | head -n 1 > {output} && printf '\n' >> {output}")
         # Add the result rows
         shell("cat {input.quast} >> {output}" 2> {log})

rule samtools_flagstat:
    input:
        expand("scratch/coverm/bamfiles/all.merged.contigs.fasta.{sample}_R1.fastq.bam", sample=config["data"])
    output:
        "results/stats/flagstat/flagstat.txt"
    conda:
        "../../../envs/samtools.yaml"
    log:
        "logs/stats/flagstat.log"
    shell:
        "samtools flagstat {input} > {output} 2> {log}"

"""
rule flagstat_convert:
    input:
        "scratch/stats/{assembler}/{treatment}/{kmers}/flagstat.txt"
    output:
        "scratch/stats/{assembler}/{treatment}/{kmers}/flagstat.linear.txt"
    params:
       run="{assembler}-{treatment}-{kmers}"
    # Create a linearized output
    # First add the assembly name as first column
    # Get the first column of the flagstat output and use tabs in stead of newlines as delimiter
    shell: "printf '{params.run}\t' > {output} && cut -d' ' -f 1 {input} | paste -s -d '\t' >> {output}"

rule flagstat_merge:
    input:
        expand("scratch/stats/{assembler}/{treatment}/{kmers}/flagstat.linear.txt", assembler=config["assembler"], treatment=config["treatment"], kmers=config["assembly-klist"]),
    output:
        "scratch/stats/flagstat.report.txt"
    run:
         # Add a header
         shell("echo 'Assembly\ttotal_reads\tsecondary\tsupplementary\tduplicates\tmapped\tpaired\tread1\tread2\tproperly_paired\twith_itself_and_mate_mapped\tsingeltons\twith_mate_mapped_different_chr\twith_mate_mapped_different_chr_q5' > {output}")
         # Add the result rows
         shell("cat {input} >> {output}")
"""
