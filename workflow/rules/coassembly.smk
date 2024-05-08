"""The rules for the coassembly part of the pipeline."""
import os

if config['sample_pooling']=='supervised':
    def get_forward_files(wildcards):
        return ["results/02_filtered_reads/" + sample + "_filt_1.fastq.gz" for sample in
                samples[samples["sample_pool"] == wildcards.sample_pool]["sample"].values]


    def get_rev_files(wildcards):
        return ["results/02_filtered_reads/" + sample + "_filt_2.fastq.gz" for sample in
                samples[samples["sample_pool"] == wildcards.sample_pool]["sample"].values]


    rule supervised_pooling:
        input:
            forward=get_forward_files,
            rev=get_rev_files
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

if config['assembly_method']=='coassembly':
    rule megahit:
        input:
            forward = "results/03_assembly/coassembly/pools/{sample_pool}_forward.fastq",
            rev = "results/03_assembly/coassembly/pools/{sample_pool}_rev.fastq",
    #        unpaired = "scratch/treatment/{treatment}_unpaired.fastq"
        output:
            contigs="results/03_assembly/coassembly/megahit_assembly/assembly_{sample_pool}/{sample_pool}_final.contigs.fa",
            # This file contains all the settings of a run. When this file is not present megahit with run in normal mode, otherwise it continues with previous settings
            # opts=protected("results/03_coassembly/{sample_pool}/{kmers}/options.json")
        # params:
        #     dir=lambda wildcards, output: os.path.dirname(output[0]),
        #     kmers = lambda wildcards: config["assembly-klist"][wildcards.kmers]
        # log: "logs/assembly/megahit/{sample_pool}/megahit.log"
        threads: 32
        conda:
            "../envs/megahit.yaml"
        shell:"megahit --out-dir results/03_assembly/coassembly/megahit_assembly/assembly_{sample_pool}/{sample_pool}_final.contigs.fa -m 0.9 -t {threads} --presets meta-large -1 {input.forward} -2 {input.rev}"

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
