rule prodigal:
    input:
        contigs=branch(config['assembly_method'] == "coassembly",
            then=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta",
            otherwise=f"{outdir}/results/03_assembly/single_sample_assembly/{{sample_pool}}/{{sample_pool}}_assembly.fasta")
    output:
        gff = f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_genes.gff",
        faa = f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_proteins.faa",
        fna = f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_orfs.fna"
    conda:
        "../envs/prodigal.yaml"
    shell:
        "prodigal -i {input.contigs} -o {output.gff} -a {output.faa} -d {output.fna} -p meta -f gff"


rule CAT:
    input:
        contigs = {rules.rename_megahit.output.fasta},
        proteins = {rules.prodigal.output.faa}
    output:
        CAT = f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/{{sample_pool}}.contig2classification.txt",
        names = f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/{{sample_pool}}.contig2classification.names.txt",
        summary = f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/{{sample_pool}}.contig2classification.names.summarise.txt"
    params:
        db = config["CAT_database"],
        tax = config["CAT_taxonomy"],
        out_dir = f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/",
        prefix = f"{{sample_pool}}",
        threads = config["threads"]
    conda:
        "../envs/cat.yaml"
    shell:
        """
        mkdir -p {params.out_dir}
        CAT_pack contigs -c {input.contigs} -n {params.threads} -d {params.db} -t {params.tax} -p {input.proteins} -o {params.out_dir}{params.prefix}
        CAT_pack add_names -i {output.CAT} -o {output.names} -t {params.tax} --only_official --exclude_scores
        CAT_pack summarise -c {input.contigs} -i {output.names} -o {output.summary}
        """

rule MetaPhlAn4:
    input:
        forward = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_1.fastq.gz",
        rev = f"{outdir}/results/02_filtered_reads/{{sample}}_filt_2.fastq.gz",
    output:
        file = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/temp_MetaPhlAn/{{sample}}.txt",
    params:
        bowtie = lambda wildcards: f"{wildcards.sample}.bowtie2.bz2",
        threads = config["threads"],
        mtphln_outdir = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/"
    conda:
        "../envs/metaphlan4.yaml"
    shell:
        """
        mkdir -p {params.mtphln_outdir}
        metaphlan {input.forward},{input.rev} --bowtie2out {params.mtphln_outdir}/{params.bowtie} --nproc {params.threads} --input_type fastq -o {output.file}
        """

rule MetaPhlAn_secondary:
    input:
        lambda wildcards: expand(f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/temp_MetaPhlAn/{{sample}}.txt", sample=samples["sample"])
    output:
        merged_table = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/merged_abundance_table.txt"
    params:
        mtphln_dir = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/",
        scripts_dir= "./workflow/scripts/"
    conda:
        "../envs/metaphlan4.yaml"
    shell:
        """
        merge_metaphlan_tables.py {input} > {output.merged_table}
        """
        
       # Rscript {params.scripts_dir}/MetaPhlAn_calculate_diversity.R -f {output.merged_table} -o {params.mtphln_dir}/beta_diversity
       # Rscript {params.scripts_dir}/MetaPhlAn_calculate_diversity.R -f {output.merged_table} -d alpha -m shannon -o {params.mtphln_dir}/alpha_diversity
       # """

rule eggnog:
    input:
        proteins = rules.prodigal.output.faa
    params:
        db = config["emapper_database"],
        threads = config["threads"],
        out_dir = f"{outdir}/results/05_prokaryote_annotation/eggnog/{{sample_pool}}/{{sample_pool}}"
    conda:
        "../envs/eggnog.yaml"
    output:
        f"{outdir}/results/05_prokaryote_annotation/eggnog/{{sample_pool}}/{{sample_pool}}.emapper.annotations"
    shell:
        "emapper.py -i {input.proteins} --cpu {params.threads} -o {params.out_dir} --data_dir {params.db} --pident 30 --query_cover 50 --subject_cover 50 --report_orthologs --override"
