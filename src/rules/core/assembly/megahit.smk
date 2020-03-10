rule megahit:
    input:
        forward = "scratch/treatment/{treatment}_forward.fastq",
        reverse = "scratch/treatment/{treatment}_reverse.fastq",
#        unpaired = "scratch/treatment/{treatment}_unpaired.fastq"
    output:
        contigs="scratch/assembly/megahit/{treatment}/{kmers}/final.contigs.fa",
        # This file contains all the settings of a run. When this file is not present megahit with run in normal mode, otherwise it continues with previous settings
        opts=protected("scratch/assembly/megahit/{treatment}/{kmers}/options.json")
    params:
        dir="scratch/assembly/megahit/{treatment}/{kmers}/",
        kmers = lambda wildcards: config["assembly-klist"][wildcards.kmers]
    log: "scratch/assembly/megahit/{treatment}/{kmers}/megahit.log"
    threads: 32
    conda:
        "../../../envs/megahit.yaml"
    shell:"megahit --continue --force --out-dir {params.dir} --tmp-dir /scratch/tmp -m 0.9 --max-read-len 302 --cpu-only -t {threads} --k-list {params.kmers} -1 {input.forward} -2 {input.reverse} 2> {log}"

rule rename_megahit:
    input:
        "scratch/assembly/megahit/{treatment}/{kmers}/final.contigs.fa"
    output:
        gzip="scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa.gz",
        fasta=temp("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa")
    run:
        shell("cat {input} | awk '{{print $1}}' | sed 's/_/contig/' > {output.fasta}")
        shell("gzip -c {output.fasta} > {output.gzip}")

