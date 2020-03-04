rule spades:
    input:
        forward = "scratch/treatment/{treatment}_forward.fastq",
        reverse = "scratch/treatment/{treatment}_reverse.fastq",
#        unpaired = "scratch/treatment/{treatment}_unpaired.fastq"
    output:
        temp("scratch/assembly/spades/{treatment}/{kmers}/contigs.fasta")
    params:
        outdir="scratch/assembly/spades/{treatment}/{kmers}/",
        kmers = lambda wildcards: config["assembly-klist"][wildcards.kmers]
    log:
        "scratch/assembly/spades/{treatment}/{kmers}/spades.log"
    threads: 32
    conda:
        "envs/spades.yaml"
    shell: "metaspades.py -m 1200 -1 {input.forward} -2 {input.reverse} --only-assembler -k {params.kmers} -t {threads} -o {params.outdir} --tmp-dir {params.outdir}/tmp/ 2>&1 > /dev/null"

rule rename_spades:
    input:
        "scratch/assembly/spades/{treatment}/{kmers}/contigs.fasta"
    output:
        gzip=protected("scratch/assembly/spades/{treatment}/{kmers}/assembly.fa.gz"),
        fasta=temp("scratch/assembly/spades/{treatment}/{kmers}/assembly.fa")
    run:
        shell("cat {input} | awk '{{print $1}}' | sed 's/_/contig/' > {output.fasta}")
        shell("gzip -c {output.fasta} > {output.gzip}")

