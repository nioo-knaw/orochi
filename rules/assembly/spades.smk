rule spades:
    input:
        forward = "{project}/treatment/{treatment}_forward.fastq",
        reverse = "{project}/treatment/{treatment}_reverse.fastq",
        unpaired = "{project}/treatment/{treatment}_unpaired.fastq"
    output:
        temp("{project}/assembly/spades/{treatment}/{kmers}/contigs.fasta")
    params:
        outdir="{project}/assembly/spades/{treatment}/{kmers}/",
        kmers = lambda wildcards: config["assembly-klist"][wildcards.kmers]
    log:
        "{project}/assembly/spades/{treatment}/{kmers}/spades.log"
    threads: 32
    conda:
        "envs/spades.yaml"
    shell: "metaspades.py -m 1200 -1 {input.forward} -2 {input.reverse} -s {input.unpaired} --only-assembler -k {params.kmers} -t {threads} -o {params.outdir} --tmp-dir {params.outdir}/tmp/ 2>&1 > /dev/null"

rule rename_spades:
    input:
        "{project}/assembly/spades/{treatment}/{kmers}/contigs.fasta"
    output:
        gzip=protected("{project}/assembly/spades/{treatment}/{kmers}/assembly.fa.gz"),
        fasta=temp("{project}/assembly/spades/{treatment}/{kmers}/assembly.fa")
    run:
        shell("cat {input} | awk '{{print $1}}' | sed 's/_/contig/' > {output.fasta}")
        shell("gzip -c {output.fasta} > {output.gzip}")

