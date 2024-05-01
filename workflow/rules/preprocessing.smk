rule fastp:
        input:
            fq1 = lambda wildcards: samples.loc[samples["sample"] == wildcards.sample].fq1.item(),
            fq2 = lambda wildcards: samples.loc[samples["sample"] == wildcards.sample].fq2.item(),

            # readF = "raw/{sample}_R1.fastq.gz",
            # readR = "raw/{sample}_R2.fastq.gz"
        output:
            cleanF = temp("results/01_trimmed_reads/{sample}_trim_1.fastq.gz"),
            cleanR = temp("results/01_trimmed_reads/{sample}_trim_2.fastq.gz"),
            report = "results/01_trimmed_reads/quality_reports/{sample}.html"
        params:
            report_name = lambda wildcards:"{wildcard.sample}"
        conda:
            "../envs/preprocessing.yaml"
        shell:
            "time fastp -i {input.fq1} -I {input.fq2} -o {output.cleanF} -O {output.cleanR} \
                        -h {output.report} -R {params.report_name} -y -l 30 -r --cut_window_size 4 \
                        --cut_mean_quality 23 --n_base_limit 0"

rule concat_host_phix:
        input:
            host = config["host_genome"],
            phix = "resources/contaminants_refs/GCF_000819615.1_ViralProj14015_genomic.fna"
        output:
            concat = temp("resources/contaminants_refs/contaminants_concat.fna")
        shell:
            "cat {input.host} {input.phix} > {output.concat}"

rule filter_host:
        input:
            readF = rules.fastp.output.cleanF,
            readR = rules.fastp.output.cleanR,
            concat = rules.concat_host_phix.output.concat,

        output:
            filterF = "results/02_filtered_reads/{sample}_filt_1.fastq.gz",
            filterR = "results/02_filtered_reads/{sample}_filt_2.fastq.gz"
        params:
            threads=config['threads']
        conda:
            "../envs/preprocessing.yaml"
        shell:
            "time bbmap.sh -Xmx52g threads={params.threads} minid=0.95 maxindel=3 \
                           in1={input.readF} in2={input.readR} \
                           outu1={output.filterF} outu2={output.filterR} \
                           ref={input.concat}"
