rule megahit:
    input:
        forward = "scratch/treatment/{treatment}_forward.fastq",
        reverse = "scratch/treatment/{treatment}_reverse.fastq",
#        unpaired = "scratch/treatment/{treatment}_unpaired.fastq"
    output:
        contigs="scratch/assembly/megahit/{treatment}/{kmers}/final.contigs.fa",
        # This file contains all the settings of a run. When this file is not present megahit with run in normal mode, otherwise it continues with previous settings
        opts=protected("scratch/assembly/megahit/{treatment}/{kmers}/opts.txt")
    params:
        dir="scratch/assembly/megahit/{treatment}/{kmers}/",
        kmers = lambda wildcards: config["assembly-klist"][wildcards.kmers]
    log: "scratch/assembly/megahit/{treatment}/{kmers}/megahit.log"
    threads: 32
    conda:
        "../../envs/megahit.yaml"
    # Parameter settings
    # meta            '--min-count 2 --k-list 21,41,61,81,99'             (generic metagenomes, default)
    # meta-sensitive  '--min-count 2 --k-list 21,31,41,51,61,71,81,91,99' (more sensitive but slower)
    # meta-large      '--min-count 2 --k-list 27,37,47,57,67,77,87'       (large & complex metagenomes, like soil)
    # bulk            '--min-count 3 --k-list 31,51,71,91,99 --no-mercy'  (experimental, standard bulk sequencing with >= 30x depth)
    # single-cell     '--min-count 3 --k-list 21,33,55,77,99,121 --merge_level 20,0.96' (experimental, single cell data)
    shell:"ulimit -m 700000000; megahit --continue --out-dir {params.dir} --tmp-dir /scratch/tmp -m 0.9 --max-read-len 302 --cpu-only -t {threads} --k-list {params.kmers} -1 {input.forward} -2 {input.reverse} 2> {log}"

rule rename_megahit:
    input:
        "scratch/assembly/megahit/{treatment}/{kmers}/final.contigs.fa"
    output:
        gzip="scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa.gz",
        fasta=temp("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa")
    run:
        shell("cat {input} | awk '{{print $1}}' | sed 's/_/contig/' > {output.fasta}")
        shell("gzip -c {output.fasta} > {output.gzip}")

