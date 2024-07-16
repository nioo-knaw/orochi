"""The rules for the coassembly part of the pipeline."""
import os

rule supervised_pooling:
    input:
        forward=lambda wildcards: get_forward_files(wildcards, "_filt_"),
        rev=lambda wildcards: get_rev_files(wildcards, "_filt_")
    output:
        forward="results/03_assembly/coassembly/pools/{sample_pool}_forward.fastq",
        rev="results/03_assembly/coassembly/pools/{sample_pool}_rev.fastq"
    run:
        shell("cat {input.forward} > {output.forward}")
        shell("cat {input.rev} > {output.rev}")

# if config['sample_pooling'] == 'simka':
#     rule simka_pooling:
#         input:
#             all_reads=lambda wildcards: expand("results/02_filtered_reads/{sample}_filt_1.fastq.gz", sample=samples[samples["sample"] == wildcards.sample]["sample"].values) + expand("results/02_filtered_reads/{sample}_filt_2.fastq.gz", sample=samples[samples["sample"] == wildcards.sample]["sample"].values)
#         output:
#             sample_pools= "results/03_coassembly/pools/{sample_pool}.fastq"

rule normal_reads:
    input:
        r1 = rules.supervised_pooling.output.forward,
        r2 = rules.supervised_pooling.output.rev
    output:
        out1="results/03_assembly/coassembly/pools/{sample_pool}_normalized_f.fastq",
        out2="results/03_assembly/coassembly/pools/{sample_pool}_normalized_r.fastq",
        hist="results/03_assembly/coassembly/pools/{sample_pool}.hist"
    params:
        kmerdepth=config['bbmap_D'],
        threads=config['threads']
    conda:
        "../envs/coassembly.yaml"
    shell:
        "bbnorm.sh target={params.kmerdepth} minprob=0.6 prefiltersize=0.50 prefilter=True min=2 in={input.r1} in2={input.r2} threads={params.threads} out={output.out1} out2={output.out2} hist={output.hist} -Xmx100g"


rule megahit:
    input:
        fwd=branch(config['normalize_reads'] == "Yes",
            then="results/03_assembly/coassembly/pools/{sample_pool}_normalized_f.fastq",
            otherwise="results/03_assembly/coassembly/pools/{sample_pool}_forward.fastq"),
        rev=branch(config['normalize_reads'] == "Yes",
            then="results/03_assembly/coassembly/pools/{sample_pool}_normalized_r.fastq",
            otherwise="results/03_assembly/coassembly/pools/{sample_pool}_rev.fastq")

    output:
        contigs="results/03_assembly/coassembly/assembly_{sample_pool}/{sample_pool}_final.contigs.fa.gz",
        # This file contains all the settings of a run. When this file is not present megahit with run in normal mode, otherwise it continues with previous settings
        # opts=protected("results/03_coassembly/{sample_pool}/{kmers}/options.json")
    params:
        # dir=lambda wildcards, output: os.path.dirname(output[0]),
        kmers = config["kmers"]
    # log: "logs/assembly/megahit/{sample_pool}/megahit.log"
    threads: 32
    conda:
        "../envs/megahit.yaml"
    shell:"megahit --out-dir results/03_assembly/coassembly/assembly_{wildcards.sample_pool} --out-prefix {wildcards.sample_pool}_final -m 0.9 --k-list {params.kmers} -t {threads} --presets meta-large -1 {input.fwd} -2 {input.rev}"

    # rule rename_megahit:
    #     input:
    #         "scratch/assembly/megahit/{sample_pool}/{kmers}/final.contigs.fa"
    #     output:
    #         gzip=temp("scratch/assembly/megahit/{sample_pool}/{kmers}/assembly.fa.gz"),
    #         fasta=temp("scratch/assembly/megahit/{sample_pool}/{kmers}/assembly.fa")
    #     log:
    #         "logs/assembly/megahit/{sample_pool}/{kmers}/rename_megahit.log"
    #     run:
    #         shell("cat {input} | awk '{{print $1}}' | sed 's/_/contig/' > {output.fasta} 2> {log}")
    #         shell("gzip -c {output.fasta} > {output.gzip}")
