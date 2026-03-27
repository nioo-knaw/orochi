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
    resources: mem_mb = 100000  # Set a high memory limit for Prodigal (100GB), but not max_mb, to still allow for parallelization
    shell:
        "prodigal -i {input.contigs} -o {output.gff} -a {output.faa} -d {output.fna} -p meta -f gff"


rule CAT:
    input:
#        contigs = {rules.rename_megahit.output.fasta},
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
    shell:
        """
        mkdir -p {params.out_dir}
        CAT_pack contigs -c {input.contigs} -n {threads} -d {params.db} -t {params.tax} -p {input.proteins} -o {params.out_dir}{params.prefix}
        CAT_pack add_names -i {output.CAT} -o {output.names} -t {params.tax} --only_official --exclude_scores
        CAT_pack summarise -c {input.contigs} -i {output.names} -o {output.summary}
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
    shell:
        """
        mkdir -p {params.mtphln_outdir}
        metaphlan {input.forward},{input.rev} --mapout {params.mtphln_outdir}/{params.bowtie} --nproc {threads} --input_type fastq -o {output.file}
        """

rule MetaPhlAn_sgb_to_gtdb:
    input:
        raw = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/temp_MetaPhlAn/{{sample}}.raw.txt"
    output:
        final = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/temp_MetaPhlAn/{{sample}}.txt"
    conda:
        "../envs/metaphlan4.yaml"
    shell:
        """
        if [ "{config[taxonomy_type]}" = "GTDB" ]; then
            sgb_to_gtdb_profile.py -i {input.raw} -o {output.final}
            rm -f {input.raw}
        else
            mv {input.raw} {output.final}
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
    shell:
        """
        merge_metaphlan_tables.py {input} {params.gtdb_flag} > {params.raw_table}
        if [ "{params.use_gtdb_reform}" = "True" ]; then
            python workflow/scripts/gtdb_reform.py \\
                -i {params.raw_table} \\
                -o {output.merged_table}
        else
            mv {params.raw_table} {output.merged_table}
        fi
        """
        
       # Rscript {params.scripts_dir}/MetaPhlAn_calculate_diversity.R -f {output.merged_table} -o {params.mtphln_dir}/beta_diversity
       # Rscript {params.scripts_dir}/MetaPhlAn_calculate_diversity.R -f {output.merged_table} -d alpha -m shannon -o {params.mtphln_dir}/alpha_diversity
       # """

rule eggnog:
    input:
        proteins = rules.prodigal.output.faa
    params:
        db = config["emapper_database"],
        out_dir = f"{outdir}/results/05_prokaryote_annotation/eggnog/{{sample_pool}}/{{sample_pool}}"
    threads:
        config["threads"]
    conda:
        "../envs/eggnog.yaml"
    output:
        raw = f"{outdir}/results/05_prokaryote_annotation/eggnog/{{sample_pool}}/{{sample_pool}}.emapper.annotations",
        adj = f"{outdir}/results/05_prokaryote_annotation/eggnog/{{sample_pool}}/{{sample_pool}}.emapper.annotations.adjusted"
    resources:
        mem_mb = 500000  # Set a high memory limit for eggNOG (500GB), but not max_mb, to still allow for parallelization
    shell:
        """
        emapper.py -i {input.proteins} --cpu {threads} -o {params.out_dir} --data_dir {params.db} --pident 30 --query_cover 50 --subject_cover 50 --report_orthologs --override
        head -n -3  <(tail -n +5 {output.raw}) > {output.adj}
        """
    # @Todo: Perhaps specify the temp dir for eggnog to avoid issues with large files?

# import pandas as pd

# df = pd.read_csv(config["samples"], sep="\t")

# SAMPLES = df["sample"].tolist()
# POOLS = sorted(df["sample_pool"].unique())

# sample_to_pool = dict(zip(df["sample"], df["sample_pool"]))
# sample_to_fq1 = dict(zip(df["sample"], df["fq1"]))
# sample_to_fq2 = dict(zip(df["sample"], df["fq2"]))

# rule salmon_assemblies1:
#     input:
#         orfs = f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_orfs.fna"
#     output:
#         index_file = f"{outdir}/results/05_prokaryote_annotation/salmon/{{sample_pool}}/{{sample_pool}}_orfs.index"
#     threads:
#         config["threads"]
#     conda:
#         "../envs/salmon.yaml"
#     resources:
#         mem_mb = 200000
#     shell:
#         """
#         salmon index -t {input.orfs} -i {output.index_file} -k 31
#         """

# rule salmon_samples2:
#     input:
#         index = lambda wc: f"salmon_index/{sample_to_pool[wc.sample]}",
#         forward = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_1.fastq.gz",
#         reverse = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_2.fastq.gz"
#     output:
#         directory(f"{outdir}/results/05_prokaryote_annotation/salmon/{{sample}}/")
#     threads:
#         config["threads"]
#     conda:
#         "../envs/salmon.yaml"
#     resources:
#         mem_mb = 200000
#     shell:
#         """
#         salmon quant -i {input.index} --libType IU -1 {input.forward} -2 {input.reverse} -p {threads} -o {output} --meta
#         """

# rule salmon_final3:
#     input:
#         quants=lambda wc: expand(
#             f"{outdir}/results/05_prokaryote_annotation/salmon/{{sample}}/",
#             sample=[s for s in SAMPLES if sample_to_pool[s] == wc.sample_pool])
#     output:
#         quant = f"{outdir}/results/05_prokaryote_annotation/salmon/{{sample_pool}}/{{sample_pool}}_ORF_TPM.tsv"
#     params:
#         sample_name = lambda wc: ",".join(
#             [s for s in SAMPLES if sample_to_pool[s] == wc.sample_pool]
#         )
#     threads:
#         config["threads"]
#     conda:
#         "../envs/salmon.yaml"
#     resources:
#         mem_mb = 200000
#     shell:
#         """
#         salmon quantmerge --quants {params.sample_name}/ --names {{params.sample_name}} --column TPM -o {output.quant}
#         """
