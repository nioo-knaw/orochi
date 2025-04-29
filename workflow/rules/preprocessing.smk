rule fastp:
        input:
            fq1 = lambda wildcards: samples.loc[samples["sample"] == wildcards.sample].fq1.item(),
            fq2 = lambda wildcards: samples.loc[samples["sample"] == wildcards.sample].fq2.item(),

            # readF = "raw/{sample}_R1.fastq.gz",
            # readR = "raw/{sample}_R2.fastq.gz"
        output:
            cleanF = temp(f"{outdir}/results/01_trimmed_reads/{{sample}}_trim_1.fastq.gz"),
            cleanR = temp(f"{outdir}/results/01_trimmed_reads/{{sample}}_trim_2.fastq.gz"),
            report_html = f"{outdir}/results/01_trimmed_reads/quality_reports/{{sample}}.html",
            report_json = f"{outdir}/results/01_trimmed_reads/quality_reports/{{sample}}.json"
        params:
            report_name = lambda wildcards:"{wildcards.sample}"
        threads:
            16 #is max nr of threads for fastp
        conda:
            "../envs/preprocessing.yaml"
        log: f"{outdir}/logs/fastp_{{sample}}.log",
        shell:
            "fastp -i {input.fq1} -I {input.fq2} -o {output.cleanF} -O {output.cleanR} \
                        -h {output.report_html} -j {output.report_json} -R {params.report_name} -w {threads} -y -l 30 -r --trim_poly_g --n_base_limit 0 2> {log}"

rule concat_host_phix:
        input:
            host = config["host_genome"],
            phix = "resources/contaminants_refs/GCF_000819615.1_ViralProj14015_genomic.fna"
        output:
            concat = temp(f"{outdir}/results/00_misc/contaminants_refs/contaminants_concat.fna")
        log: f"{outdir}/logs/concat_host_phix.log"
        shell:
            "cat {input.host} {input.phix} > {output.concat} 2> {log}"

rule build_index:
    conda:
        "../envs/preprocessing.yaml"
    input: 
        reference=f"{outdir}/results/00_misc/contaminants_refs/contaminants_concat.fna"
    output: 
        ref_index=directory(f"{outdir}/results/00_misc/contaminants_refs/ref/")
    params:
        # threads=config['threads'],
        memory=config['bbmap_mem']
    resources:
        mem_mb=1500000
    threads:
        workflow.cores * 0.5
    log: f"{outdir}/logs/build_index.log"
    shell:
        "bbmap.sh ref={input.reference} path={output.ref_index} threads={threads} {params.memory} 2> {log}"

rule filter_host:
        input:
            readF = rules.fastp.output.cleanF,
            readR = rules.fastp.output.cleanR,
            concat = f"{outdir}/results/00_misc/contaminants_refs/contaminants_concat.fna",
            ref_index=f"{outdir}/results/00_misc/contaminants_refs/ref/"

        output:
            filterF = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_1.fastq.gz",
            filterR = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_2.fastq.gz"
        params:
            # threads=config['threads'],
            memory=config['bbmap_mem'],
            ref_dir=f"{outdir}/results/00_misc/contaminants_refs"
        threads:
            workflow.cores * 0.5
        resources:
            mem_mb=1500000
        conda:
            "../envs/preprocessing.yaml"
        log: f"{outdir}/logs/filter_host_{{sample}}.log"
        shell:
            "bbmap.sh threads={threads} minid=0.95 maxindel=3 \
                           in1={input.readF} in2={input.readR} path={input.ref_index} \
                           outu1={output.filterF} outu2={output.filterR} {params.memory} 2> {log}"



