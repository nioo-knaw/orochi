
rule spades:
    input:
        forward = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_1.fastq.gz",
        rev = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_2.fastq.gz",
    output:
        f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}/contigs.fasta",
    params:
        outdir = f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}",
        kmers = config['kmers']
    threads: int(workflow.cores * 0.8)
    resources:
        mem_mb = config['max_mem']
    conda:
        "../envs/single_assembly.yaml"
    log:
        f"{outdir}/logs/spades/spades_{{sample}}.log"
    shell:
        """
        spades.py \
            --meta \
            -m 1200 \
            -1 {input.forward} \
            -2 {input.rev} \
            --only-assembler \
            -k {params.kmers} \
            -t {threads} \
            -o {params.outdir} \
            --tmp-dir {params.outdir}/tmp/spades \
            > {log} 2>&1
        """


rule rename_spades:
    input:
        contigs = rules.spades.output
    output:
        gzip = f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}/{{sample}}_assembly.fasta.gz",
        fasta = temp(f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}/{{sample}}_assembly.fasta")
    run:
        shell("cat {input.contigs} | awk '{{print $1}}' | sed 's/NODE/contig/' > {output.fasta}")
        shell("gzip -c {output.fasta} > {output.gzip}")

rule assembly_quality_single:
    input:
        assembly = rules.rename_spades.output.gzip
    output:
        mq_out = f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}/quast_results/report.html"
    params:
        # threads = config['threads'],
        outdir = f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}/quast_results/"
    threads:
        int(workflow.cores * 0.8)
    conda:
        "../envs/single_assembly.yaml"
    log:
        f"{outdir}/logs/metaquast/metaquast_{{sample}}.log"
    shell:
        """
        metaquast.py \
        {input.assembly} \
        --no-krona \
        --threads {threads} \
        -o {params.outdir}\
        > {log} 2>&1
        """

# Note: single-sample binning coverage now lives in binning.smk as
# coverm_coverage (maps against the size-filtered contigs binning actually
# uses, and writes MetaBAT-compatible depth output directly).
