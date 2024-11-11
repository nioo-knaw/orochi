""" Rules related to binning/creating MAGs """

import pandas as pd
sample_to_pool = pd.read_csv(config['samples'], sep='\t', index_col="sample")["sample_pool"].to_dict()

rule index_bwa:
    input:
        assembly=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta"
    output:
        f"{outdir}/results/06_binning/mapping/{{sample_pool}}/{{sample_pool}}_assembly.fasta.bwt"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input.assembly}"

rule mapping:
    input:
        assembly=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta",
        assembly_index=f"{outdir}/results/06_binning/mapping/{{sample_pool}}/{{sample_pool}}_assembly.fasta.bwt",
        forward=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_normalized_f.fastq",
        rev=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_normalized_r.fastq"

    output:
        bam=f"{outdir}/results/06_binning/mapping/{{sample_pool}}/{{sample_pool}}_sorted.bam",
        bam_index=f"{outdir}/results/06_binning/mapping/{{sample_pool}}/{{sample_pool}}_sorted.bam.bai",
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa mem -t 4 {input.assembly} {input.forward} {input.rev} | samtools sort -@ 4 -o {output.bam}"
        "samtools index {output.bam}"

# rule mapping:
#     input:
#         assembly=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta",
#         forward=lambda wildcards: get_forward_files(wildcards, "_filt_"),
#         rev=lambda wildcards: get_rev_files(wildcards, "_filt_"),
#
#     output:
#         bam=f"{outdir}/results/06_binning/mapping/{{sample_pool}}/{{sample_pool}}_{{sample}}_sorted.bam",
#         # bam_index=f"{outdir}/results/06_binning/mapping/{{sample_pool}}/{{sample_pool}}_{{sample}}_sorted.bam.bai",
#     conda:
#         "../envs/bwa.yaml"
#     shell:
#         """
#         for i in $(seq 0 $(expr {{input.forward|length}} - 1)); do
#             forward=${{input.forward[i]}}
#             reverse=${{input.rev[i]}}
#             echo "Processing $forward and $reverse"
#             bwa mem -t 4 {input.assembly} "$forward" "$reverse" | \
#                 samtools sort -@ 4 -o {output.bam}
#         done
#         """
        # "bwa mem -t 4 {input.assembly} {input.forward} {input.rev} | samtools sort -@ 4 -o {output.bam}"
        # "samtools index {output.bam}"

rule jgi_summarize_bam_contig_depths:
    input:
        f"{outdir}/results/06_binning/mapping/{{sample}}/{{sample}}_sorted.bam"
    output:
        f"{outdir}/results/06_binning/mapping/{{sample}}/{{sample}}_depth.txt"
    conda:
        "../envs/metabat2.yaml"
    shell:
        "jgi_summarize_bam_contig_depths --outputDepth {output} {input}"

rule metabat2:
    input:
        sorted_bam=f"{outdir}/results/06_binning/mapping/{{sample}}/{{sample}}_sorted.bam",
        assembly=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample}}/{{sample}}_assembly.fasta",
        depth=f"{outdir}/results/06_binning/mapping/{{sample}}/{{sample}}_depth.txt"
    output:
        bin_dir=directory(f"{outdir}/results/06_binning/metabat2/{{sample}}/{{sample}}_bins/"),
        bin_tsv=f"{outdir}/results/06_binning/metabat2/{{sample}}/{{sample}}_bins.tsv"

    conda:
        "../envs/metabat2.yaml"
    shell:
        "metabat2 -i {input.assembly} -a {input.depth} -o {output.bin_dir} -t 4"

rule maxbin2:
    input:
        f"{outdir}/results/06_binning/{{sample}}/{{sample}}_sorted.bam",
        f"{outdir}/results/03_assembly/coassembly/assembly_{{sample}}/{{sample}}_assembly.fasta",
    output:
        bins=f"{outdir}/results/06_binning/{{sample}}/{{sample}}_bins"
    conda:
        "../envs/maxbin2.yaml"
    shell:
        "run_MaxBin.pl -contig {input[1]} -reads {input[0]} -out {output} -thread 4"

rule vamb:
    input:
        f"{outdir}/results/06_binning/{{sample}}/{{sample}}_sorted.bam",
        f"{outdir}/results/03_assembly/coassembly/assembly_{{sample}}/{{sample}}_assembly.fasta",
    output:
        bins=f"{outdir}/results/06_binning/{{sample}}/{{sample}}_bins"
    conda:
        "../envs/vamb.yaml"
    shell:
        "vamb --outdir {output} --fasta {input[1]} --bam {input[0]} --minfasta 2000 --minreads 2 --threads 4"