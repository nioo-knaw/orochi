if config['assembler']=='megahit':
    rule megahit:
        input:
            forward = "scratch/treatment/{treatment}_forward.fastq",
            rev = "scratch/treatment/{treatment}_rev.fastq",
    #        unpaired = "scratch/treatment/{treatment}_unpaired.fastq"
        output:
            contigs="scratch/assembly/megahit/{treatment}/{kmers}/final.contigs.fa",
            # This file contains all the settings of a run. When this file is not present megahit with run in normal mode, otherwise it continues with previous settings
            opts=protected("scratch/assembly/megahit/{treatment}/{kmers}/options.json")
        params:
            dir=lambda wildcards, output: os.path.dirname(output[0]),
            kmers = lambda wildcards: config["assembly-klist"][wildcards.kmers]
        log: "logs/assembly/megahit/{treatment}/{kmers}/megahit.log"
        threads: 80
        conda:
            "../../../envs/megahit.yaml"
        shell:"megahit --continue --force --out-dir {params.dir} --tmp-dir /tmp -m 0.9 -t {threads} --k-list {params.kmers} -1 {input.forward} -2 {input.rev} 2> {log}"

    rule rename_megahit:
        input:
            "scratch/assembly/megahit/{treatment}/{kmers}/final.contigs.fa"
        output:
            gzip="scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa.gz",
            fasta=temp("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa")
        log:
            "logs/assembly/megahit/{treatment}/{kmers}/rename_megahit.log"
        run:
            shell("cat {input} | awk '{{print $1}}' | sed 's/_/contig/' > {output.fasta} 2> {log}")
            shell("gzip -c {output.fasta} > {output.gzip}")

