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
        conda:
            "../envs/preprocessing.yaml"
        log: f"{outdir}/logs/fastp_{{sample}}.log",
        shell:
            "time fastp -i {input.fq1} -I {input.fq2} -o {output.cleanF} -O {output.cleanR} \
                        -h {output.report_html} -j {output.report_json} -R {params.report_name} -y -l 30 -r --trim_poly_g --n_base_limit 0 2> {log}"

rule concat_host_phix:
        input:
            host = config["host_genome"],
            phix = "resources/contaminants_refs/GCF_000819615.1_ViralProj14015_genomic.fna"
        output:
            concat = temp("resources/contaminants_refs/contaminants_concat.fna")
        shell:
            "cat {input.host} {input.phix} > {output.concat}"

rule build_index:
    conda:
        "../envs/preprocessing.yaml"
    input: 
        reference="resources/contaminants_refs/contaminants_concat.fna"
    output: 
        ref_index=directory("ref/")
    params:
        threads=config['threads'],
        memory=config['bbmap_mem']
    shell:
        "bbmap.sh ref={input.reference} threads={params.threads} {params.memory}"

rule filter_host:
        input:
            readF = rules.fastp.output.cleanF,
            readR = rules.fastp.output.cleanR,
            concat = "resources/contaminants_refs/contaminants_concat.fna",
            ref_index="ref/"

        output:
            filterF = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_1.fastq.gz",
            filterR = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_2.fastq.gz"
        params:
            threads=config['threads'],
            memory=config['bbmap_mem']
        conda:
            "../envs/preprocessing.yaml"
        shell:
            "time bbmap.sh threads={params.threads} minid=0.95 maxindel=3 \
                           in1={input.readF} in2={input.readR} \
                           outu1={output.filterF} outu2={output.filterR} {params.memory}"



