""" Rules related to binning/creating MAGs """

# Estimating coverage with fairy. First sketch the reads.
rule fairy_sketch:
    input:
        forward=f"{outdir}/results/02_filtered_reads/{{sample}}_filt_1.fastq.gz",
        rev=f"{outdir}/results/02_filtered_reads/{{sample}}_filt_2.fastq.gz"
    output:
        sketch_file=f"{outdir}/results/06_binning/coverage/fairy_sketch/{{sample}}_filt_1.fastq.gz.paired.bcsp",
    params:
        sketch_dir=f"{outdir}/results/06_binning/coverage/fairy_sketch"
    conda:
        "../envs/fairy.yaml"
    shell:
        "fairy sketch -1 {input.forward} -2 {input.rev} -d {params.sketch_dir}"

# Estimate coverage with fairy, using sketched reads and assembly.
rule fairy_coverage:
    input:
        assembly=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta",
        sketch_files=expand(f"{outdir}/results/06_binning/coverage/fairy_sketch/{{sample}}_filt_1.fastq.gz.paired.bcsp", sample=samples["sample"])
    output:
        coverage_file=f"{outdir}/results/06_binning/coverage/coverage_{{sample_pool}}.tsv"
    params:
        threads=config['threads']
    conda:
        "../envs/fairy.yaml"
    shell:
        "fairy coverage {input.sketch_files} {input.assembly} -t {params.threads} -o {output.coverage_file}"

# rule index_bwa:
#     input:
#         assembly=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta"
#     output:
#         f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta.bwt"
#     conda:
#         "../envs/bwa.yaml"
#     shell:
#         "bwa index {input.assembly}"

# rule mapping:
#     input:
#         assembly=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta",
#         assembly_index=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta.bwt",
#         forward=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_normalized_f.fastq",
#         rev=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_normalized_r.fastq"
#
#     output:
#         bam=f"{outdir}/results/06_binning/mapping/{{sample_pool}}/{{sample_pool}}_sorted.bam",
#         bam_index=f"{outdir}/results/06_binning/mapping/{{sample_pool}}/{{sample_pool}}_sorted.bam.bai",
#     params:
#         threads=config['threads']
#     conda:
#         "../envs/bwa.yaml"
#     shell:
#         "bwa mem -t {params.threads} {input.assembly} {input.forward} {input.rev} | samtools sort -@ {params.threads} -o {output.bam}"
#         "samtools index {output.bam}"
#
#
# rule jgi_summarize_bam_contig_depths:
#     input:
#         f"{outdir}/results/06_binning/mapping/{{sample}}/{{sample}}_sorted.bam"
#     output:
#         f"{outdir}/results/06_binning/mapping/{{sample}}/{{sample}}_depth.txt"
#     conda:
#         "../envs/metabat2.yaml"
#     shell:
#         "jgi_summarize_bam_contig_depths --outputDepth {output} {input}"

rule metabat2:
    input:
        assembly=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta",
        depth=f"{outdir}/results/06_binning/coverage/coverage_{{sample_pool}}.tsv"
    output:
        bin_dir=directory(f"{outdir}/results/06_binning/metabat2/{{sample_pool}}/{{sample_pool}}_bins/{{sample_pool}}_metabat2"),
        completed=f"{outdir}/results/06_binning/metabat2/{{sample_pool}}/{{sample_pool}}_metabat2.done"
    params:
        threads=config['threads']

    conda:
        "../envs/metabat2.yaml"
    shell:
        "metabat2 -i {input.assembly} -a {input.depth} -o {output.bin_dir} -t {params.threads} && touch {output.completed}"

# rule maxbin2:
#     input:
#         f"{outdir}/results/06_binning/{{sample}}/{{sample}}_sorted.bam",
#         f"{outdir}/results/03_assembly/coassembly/assembly_{{sample}}/{{sample}}_assembly.fasta",
#     output:
#         bins=f"{outdir}/results/06_binning/{{sample}}/{{sample}}_bins"
#     conda:
#         "../envs/maxbin2.yaml"
#     shell:
#         "run_MaxBin.pl -contig {input[1]} -reads {input[0]} -out {output} -thread 4"
#
# rule vamb:
#     input:
#         f"{outdir}/results/06_binning/{{sample}}/{{sample}}_sorted.bam",
#         f"{outdir}/results/03_assembly/coassembly/assembly_{{sample}}/{{sample}}_assembly.fasta",
#     output:
#         bins=f"{outdir}/results/06_binning/{{sample}}/{{sample}}_bins"
#     conda:
#         "../envs/vamb.yaml"
#     shell:
#         "vamb --outdir {output} --fasta {input[1]} --bam {input[0]} --minfasta 2000 --minreads 2 --threads 4"