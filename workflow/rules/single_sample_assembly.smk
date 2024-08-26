
rule spades:
    input:
        forward = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_1.fastq.gz",
        rev = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_2.fastq.gz",
    output:
        f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}/{{sample}}_contigs.fasta",
    params:
        outdir = f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}",
        kmers = config['kmers'],
    conda:
        "../envs/single_assembly.yaml"
    shell: "spades.py --meta -m 1200 -1 {input.forward} -2 {input.rev} --only-assembler -k {params.kmers} -t {threads} -o {params.outdir} --tmp-dir {params.outdir}/tmp/"

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
        threads = config['threads'],
        outdir = f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}/quast_results/"
    conda:
        "../envs/single_assembly.yaml"
    shell: "metaquast.py {input.assembly} --no-icarus --threads {params.threads} -o {params.outdir}"

rule coverm:
    input:
        contigs_f = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_1.fastq.gz",
        contigs_r = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_2.fastq.gz",
        assembly = rules.rename_spades.output.gzip
    output:
        coverm_out = f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}/quality/coverage.tsv"
    params:
        threads = config['threads']
    conda:
        "../envs/single_assembly.yaml"
    shell: "coverm contig --mapper bwa-mem --reference {input.assembly} -1 {input.contigs_f} -2 {input.contigs_r} --threads {params.threads} > {output.coverm_out}"
