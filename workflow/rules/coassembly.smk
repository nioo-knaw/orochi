"""The rules for the coassembly part of the pipeline."""
import os

rule supervised_pooling:
    input:
        samples=branch(config['sample_pooling'] == "unsupervised",
            then=rules.clustering.output.smp,
            otherwise="config/samples.tsv"), # rules.clustering.output.smp,

        forward=lambda wildcards: get_forward_files(wildcards, "_filt_"),
        rev=lambda wildcards: get_rev_files(wildcards, "_filt_")
    output:
        forward=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
        rev=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_rev.fastq.gz"
    run:
        shell("cat {input.forward} > {output.forward}")
        shell("cat {input.rev} > {output.rev}")


rule normal_reads:
    input:
        r1 = rules.supervised_pooling.output.forward,
        r2 = rules.supervised_pooling.output.rev
    output:
        out1=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_normalized_f.fastq",
        out2=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_normalized_r.fastq",
        hist=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}.hist"
    params:
        kmerdepth=config['bbmap_D'],
        threads=config['threads'],
        memory=config['bbmap_mem']
    benchmark:
        f"{outdir}/results/benchmark/normal_reads/{{sample_pool}}.tsv"
    conda:
        "../envs/coassembly.yaml"
    shell:
        "bbnorm.sh target={params.kmerdepth} minprob=0.6 prefiltersize=0.50 prefilter=True min=2 in={input.r1} in2={input.r2} threads={params.threads} out={output.out1} out2={output.out2} hist={output.hist} {params.memory}"


rule megahit:
    input:
        fwd=branch(config['normalize_reads'] == "Yes",
            then=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_normalized_f.fastq",
            otherwise=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq"),
        rev=branch(config['normalize_reads'] == "Yes",
            then=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_normalized_r.fastq",
            otherwise=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_rev.fastq")

    output:
        contigs=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_final.contigs.fa",
        # This file contains all the settings of a run. When this file is not present megahit with run in normal mode, otherwise it continues with previous settings
        # opts=protected(f"{outdir}/results/03_coassembly/{{sample_pool}}/{{params.kmers}}/options.json")
    params:
        # dir=lambda wildcards, output: os.path.dirname(output[0]),
        kmers = config["kmers"],
        output_dir = f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}",
        memory=config['megahit_mem']
    # log: "{outdir}/logs/assembly/megahit/{sample_pool}/megahit.log"
    benchmark:
        f"{outdir}/results/benchmark/megahit/{{sample_pool}}.tsv"
    threads: 32
    conda:
        "../envs/megahit.yaml"
    shell:"megahit -f --out-dir {params.output_dir} --out-prefix {wildcards.sample_pool}_final -m {params.memory} --k-list {params.kmers} -t {threads} --presets meta-large -1 {input.fwd} -2 {input.rev}"

rule rename_megahit:
    input:
        rules.megahit.output
        
    output:
        fasta=temp(f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta"),
        gzip=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta.gz",
        done= temp(f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.done")
    run:
        shell("cat {input} | awk '{{print $1}}' | sed 's/_/contig/' > {output.fasta}")
        shell("gzip -c {output.fasta} > {output.gzip}")
        shell("touch {output.done}")

