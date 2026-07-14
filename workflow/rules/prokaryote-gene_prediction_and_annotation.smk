rule prodigal:
    input:
         contigs=f"{outdir}/results/03_assembly/size_filtered/{{sample_pool}}_{minsize}/contigs_{{sample_pool}}_{minsize}.fasta"
#        contigs=branch(config['assembly_method'] == "coassembly",
#            then=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta",
#            otherwise=f"{outdir}/results/03_assembly/single_sample_assembly/{{sample_pool}}/{{sample_pool}}_assembly.fasta")
    output:
        gff = f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_genes.gff",
        faa = f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_proteins.faa",
        fna = f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_orfs.fna"
    conda:
        "../envs/prodigal.yaml"
    log: f"{outdir}/logs/prodigal/prodigal_{{sample_pool}}.log"
    resources: mem_mb = 100000  # Set a high memory limit for Prodigal (100GB), but not max_mb, to still allow for parallelization
    shell:
        """
        prodigal \
            -i {input.contigs} \
            -o {output.gff} \
            -a {output.faa} \
            -d {output.fna} \
            -p meta \
            -f gff \
            > {log} 2>&1
            """


rule CAT:
    input:
        contigs = f"{outdir}/results/03_assembly/size_filtered/{{sample_pool}}_{minsize}/contigs_{{sample_pool}}_{minsize}.fasta",
        proteins = {rules.prodigal.output.faa}
    output:
        CAT = f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/{{sample_pool}}.contig2classification.txt",
        names = f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/{{sample_pool}}.contig2classification.names.txt",
        summary = f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/{{sample_pool}}.contig2classification.names.summarise.txt",
        alignment = f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/{{sample_pool}}.alignment.diamond"
    params:
        db = config["CAT_database"],
        tax = config["CAT_taxonomy"],
        out_dir = f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/",
        prefix = f"{{sample_pool}}"
    threads:
        config["threads"]
    resources:
        mem_mb = 100000  # Set a high memory limit for CAT (100GB), but not max_mb, to still allow for parallelization
    conda:
        "../envs/cat.yaml"
    log:
        f"{outdir}/logs/CAT/CAT_{{sample_pool}}.log"
    shell:
        """
        mkdir -p {params.out_dir}
        CAT_pack contigs \
            -c {input.contigs} \
            -n {threads} \
            -d {params.db} \
            -t {params.tax} \
            -p {input.proteins} \
            -o {params.out_dir}{params.prefix} \
            --force \
            > {log} 2>&1
        CAT_pack add_names \
            -i {output.CAT} \
            -o {output.names} \
            -t {params.tax} \
            --only_official \
            --exclude_scores \
            --force \
            >> {log} 2>&1
        CAT_pack summarise \
            -c {input.contigs} \
            -i {output.names} \
            -o {output.summary} \
            --force \
            >> {log} 2>&1
        """

rule MetaPhlAn4:
    input:
        forward = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_1.fastq.gz",
        rev = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_2.fastq.gz",
    output:
        file = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/temp_MetaPhlAn/{{sample}}.raw.txt",
    params:
        bowtie = lambda wildcards: f"{wildcards.sample}.bowtie2.bz2",
        mtphln_outdir = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/"
    threads:
        64
    resources:
        mem_mb = 500000  # Still allows for parallelization, but sets a high memory limit for MetaPhlAn (500GB)
    conda:
        "../envs/metaphlan4.yaml"
    log:
        f"{outdir}/logs/MetaPhlAn4/metaplhan_{{sample}}.log"
    shell:
        """
        mkdir -p {params.mtphln_outdir}
        metaphlan {input.forward},{input.rev} --mapout {params.mtphln_outdir}/{params.bowtie} \
            --nproc {threads} --input_type fastq -o {output.file} \
            > {log} 2>&1
        """

rule MetaPhlAn_sgb_to_gtdb:
    input:
        raw = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/temp_MetaPhlAn/{{sample}}.raw.txt"
    output:
        final = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/temp_MetaPhlAn/{{sample}}.txt"
    conda:
        "../envs/metaphlan4.yaml"
    log:
        f"{outdir}/logs/MetaPhlAn4/MetaPhlAn_sgb_to_gtdb.log"
    shell:
        """
        if [ "{config[taxonomy_type]}" = "GTDB" ]; then
            sgb_to_gtdb_profile.py -i {input.raw} -o {output.final} 2>> {log}
            rm -f {input.raw} 2>> {log}
        else
            mv {input.raw} {output.final} 2>> {log}
        fi
        """

rule MetaPhlAn_secondary:
    input:
        expand(f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/temp_MetaPhlAn/{{sample}}.txt",
            sample = samples["sample"])
    output:
        merged_table = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/merged_abundance_table.txt"
    params:
        mtphln_dir = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/",
        gtdb_flag = lambda wc: (
            "--gtdb_profiles"
            if config["taxonomy_type"] == "GTDB"
            else ""
        ),
        raw_table      = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/merged_abundance_table.raw.txt",
        use_gtdb_reform = config["taxonomy_type"] == "GTDB",
#       scripts_dir= "./workflow/scripts/"
    conda:
        "../envs/metaphlan4.yaml"
    log:
        f"{outdir}/logs/MetaPhlAn4/MetaPhlAn_secondary.log"
    shell:
        """
        merge_metaphlan_tables.py {input} {params.gtdb_flag} > {params.raw_table} 2> {log}
        
        if [ "{params.use_gtdb_reform}" = "True" ]; then
            python workflow/scripts/gtdb_reform.py \\
                -i {params.raw_table} \\
                -o {output.merged_table} \
                2 >> {log}
        else
            mv {params.raw_table} {output.merged_table} 2>> {log}
        fi
        """
        
       # Rscript {params.scripts_dir}/MetaPhlAn_calculate_diversity.R -f {output.merged_table} -o {params.mtphln_dir}/beta_diversity
       # Rscript {params.scripts_dir}/MetaPhlAn_calculate_diversity.R -f {output.merged_table} -d alpha -m shannon -o {params.mtphln_dir}/alpha_diversity
       # """

rule eggnog:
    input:
        proteins = rules.prodigal.output.faa
    output:
        raw = f"{outdir}/results/05_prokaryote_annotation/eggnog/{{sample_pool}}/{{sample_pool}}.emapper.annotations",
        adj = f"{outdir}/results/05_prokaryote_annotation/eggnog/{{sample_pool}}/{{sample_pool}}.emapper.annotations.adjusted"
    params:
        db = config["emapper_database"],
        out_dir = f"{outdir}/results/05_prokaryote_annotation/eggnog/{{sample_pool}}/{{sample_pool}}",
        temp_dir = lambda wildcards: os.path.join(config["tmpdir"], "eggnog", wildcards.sample_pool)
    threads:
        config["threads"]
    conda:
        "../envs/eggnog.yaml"
    log:
        f"{outdir}/logs/eggnog/eggnog_{{sample_pool}}.log"
    resources:
        mem_mb = 500000  # Set a high memory limit for eggNOG (500GB), but not max_mb, to still allow for parallelization
    shell:
        """
        mkdir -p {params.temp_dir}
        emapper.py \
            -i {input.proteins} \
            --cpu {threads} \
            -o {params.out_dir} \
            --temp_dir {params.temp_dir} \
            --data_dir {params.db} \
            --pident 30 \
            --query_cover 50 \
            --subject_cover 50 \
            --report_orthologs \
            --override \
            > {log} 2>&1
        head -n -3  <(tail -n +5 {output.raw}) > {output.adj}
        """

import pandas as pd

df = pd.read_csv(config["samples"], sep="\t")

SAMPLES = df["sample"].tolist()
POOLS = sorted(df["sample_pool"].unique())

sample_to_pool = dict(zip(df["sample"], df["sample_pool"]))
sample_to_fq1 = dict(zip(df["sample"], df["fq1"]))
sample_to_fq2 = dict(zip(df["sample"], df["fq2"]))

ASSEMBLY_UNITS = POOLS if config["assembly_method"] == "coassembly" else SAMPLES

def assembly_unit_for_sample(sample):
    if config["assembly_method"] == "coassembly":
        return sample_to_pool[sample]
    return sample

def samples_for_assembly_unit(assembly_unit):
    if config["assembly_method"] == "coassembly":
        return [s for s in SAMPLES if sample_to_pool[s] == assembly_unit]
    return [assembly_unit]

rule salmon_assemblies1:
    input:
        orfs = f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_orfs.fna"
    output:
        index_file = directory(f"{outdir}/results/05_prokaryote_annotation/salmon/indexes/{{sample_pool}}/{{sample_pool}}_orfs.index")
    threads:
        config["threads"]
    conda:
        "../envs/salmon.yaml"
    log:
        f"{outdir}/logs/salmon_assemblies1/salmon_assemblies1_{{sample_pool}}.log"
    resources:
        mem_mb = 200000
    shell:
        """
        salmon index -t {input.orfs} -i {output.index_file} -k 31 \
            2> {log}
        """

rule salmon_samples2:
    input:
        index=lambda wc: f"{outdir}/results/05_prokaryote_annotation/salmon/indexes/{assembly_unit_for_sample(wc.sample)}/{assembly_unit_for_sample(wc.sample)}_orfs.index",
        forward=f"{outdir}/results/02_filtered_reads/{{sample}}_filt_1.fastq.gz",
        rev=f"{outdir}/results/02_filtered_reads/{{sample}}_filt_2.fastq.gz"
    output:
        directory(f"{outdir}/results/05_prokaryote_annotation/salmon/quants/{{sample}}")

    threads:
        config["threads"]
    conda:
        "../envs/salmon.yaml"
    log:
        f"{outdir}/logs/salmon_samples2/salmon_samples2_{{sample}}.log"
    resources:
        mem_mb = 200000
    shell:
        """
        salmon quant -i {input.index} --libType IU -1 {input.forward} -2 {input.rev} -p {threads} -o {output} --meta \
            2> {log}
        """

rule salmon_final3:
    input:
        quants=lambda wc: expand(
            f"{outdir}/results/05_prokaryote_annotation/salmon/quants/{{sample}}/",
            sample=samples_for_assembly_unit(wc.sample_pool))
    output:
        quant = f"{outdir}/results/05_prokaryote_annotation/salmon/merged/{{sample_pool}}/{{sample_pool}}_ORF_TPM.tsv"
    params:
        sample_dir = lambda wc: " ".join(
            [f"{outdir}/results/05_prokaryote_annotation/salmon/quants/{s}"
             for s in samples_for_assembly_unit(wc.sample_pool)]
        ),
        sample_name = lambda wc: " ".join(
            samples_for_assembly_unit(wc.sample_pool)
        ),
    threads:
        config["threads"]
    conda:
        "../envs/salmon.yaml"
    log:
        f"{outdir}/logs/salmon_final3/salmon_final3_{{sample_pool}}.log"
    resources:
        mem_mb = 200000
    shell:
        """
        salmon quantmerge --quants {params.sample_dir} --names {params.sample_name} --column TPM -o {output.quant} \
            2> {log}
        """
