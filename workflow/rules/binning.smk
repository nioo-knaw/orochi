""" Rules related to binning/creating MAGs """
import pandas as pd

# @Todo: Figure out a way to also do binning for single-sample assembly. We need to use coverM here. Maybe different rules for single-sample assembly?

binning_tools = {
    "metabat2": "fa",
    "maxbin2": "fasta"
}

# Estimating coverage with fairy. First sketch the reads.
rule fairy_sketch:
    input:
        forward=f"{outdir}/results/02_filtered_reads/{{sample}}_filt_1.fastq.gz",
        rev=f"{outdir}/results/02_filtered_reads/{{sample}}_filt_2.fastq.gz"
    output:
        sketch_file=f"{outdir}/results/06_binning/coverage/fairy_sketch/{{sample}}_filt_1.fastq.gz.paired.bcsp"
    params:
        sketch_dir=f"{outdir}/results/06_binning/coverage/fairy_sketch"
    resources:
        mem_mb=lambda wildcards: int(config['max_mem'] * 0.4)
    conda:
        "../envs/fairy.yaml"
    shell:
        "fairy sketch -1 {input.forward} -2 {input.rev} -d {params.sketch_dir}"

# Estimate coverage with fairy, using sketched reads and assembly.
rule fairy_coverage:
    input:
#        assembly=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta",
        assembly=f"{outdir}/results/03_assembly/size_filtered/{{sample_pool}}_{minsize}/contigs_{{sample_pool}}_{minsize}.fasta",
        sketch_files=expand(f"{outdir}/results/06_binning/coverage/fairy_sketch/{{sample}}_filt_1.fastq.gz.paired.bcsp", sample=samples["sample"])
    output:
        coverage_file=f"{outdir}/results/06_binning/coverage/fairy/coverage_{{sample_pool}}.tsv"
    threads:
        config['threads']
    resources:
        mem_mb=config['max_mem']
    conda:
        "../envs/fairy.yaml"
    shell:
        "fairy coverage {input.sketch_files} {input.assembly} -t {threads} -o {output.coverage_file}"


checkpoint metabat2:
    input:
#        assembly=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta",
        assembly=f"{outdir}/results/03_assembly/size_filtered/{{sample_pool}}_{minsize}/contigs_{{sample_pool}}_{minsize}.fasta",

        depth=f"{outdir}/results/06_binning/coverage/fairy/coverage_{{sample_pool}}.tsv"
    output:
        bin_dir=directory(f"{outdir}/results/06_binning/metabat2/{{sample_pool}}/{{sample_pool}}_bins")
        # completed=f"{outdir}/results/06_binning/metabat2/{{sample_pool}}/{{sample_pool}}_metabat2.done"
    params:
        bin_prefix=f"{outdir}/results/06_binning/metabat2/{{sample_pool}}/{{sample_pool}}_bins/{{sample_pool}}_bin"
    threads:
        config['threads']
    resources:
        mem_mb=config['max_mem']
    conda:
        "../envs/metabat2.yaml"
    shell:
        "metabat2 -i {input.assembly} -a {input.depth} -o {params.bin_prefix} -t {threads}"


def make_maxbin_coverage(input_file, output_file):
    """Filter columns with '-var' in their header name."""
    df = pd.read_csv(input_file, sep='\t')
    filtered_df = df.loc[:, ~df.columns.str.contains('-var')]
    filtered_df = filtered_df.loc[:, ~filtered_df.columns.str.contains('contigLen')]
    filtered_df = filtered_df.loc[:, ~filtered_df.columns.str.contains('totalAvgDepth')]
    filtered_df.to_csv(output_file, sep='\t', index=False, header=False)

rule maxbin_coverage:
    input:
        fairy_input=f"{outdir}/results/06_binning/coverage/fairy/coverage_{{sample_pool}}.tsv"
    output:
        maxbin_coverage=f"{outdir}/results/06_binning/coverage/fairy/maxbin2/coverage_{{sample_pool}}_maxbin.tsv"
    run:
        make_maxbin_coverage(input_file=input.fairy_input, output_file=output.maxbin_coverage)


checkpoint maxbin2:
    input:
#        assembly=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta",
        assembly=f"{outdir}/results/03_assembly/size_filtered/{{sample_pool}}_{minsize}/contigs_{{sample_pool}}_{minsize}.fasta",
        coverage=f"{outdir}/results/06_binning/coverage/fairy/maxbin2/coverage_{{sample_pool}}_maxbin.tsv"
    output:
        bin_dir=directory(f"{outdir}/results/06_binning/maxbin2/{{sample_pool}}/{{sample_pool}}_bins")
        # summary=f"{outdir}/results/06_binning/maxbin2/{{sample_pool}}/{{sample_pool}}_bins/{{sample_pool}}_bin.summary",
        # marker_counts=f"{outdir}/results/06_binning/maxbin2/{{sample_pool}}/{{sample_pool}}_bins/{{sample_pool}}_bin.marker",
        # marker_genes=f"{outdir}/results/06_binning/maxbin2/{{sample_pool}}/{{sample_pool}}_bins/{{sample_pool}}_bin.marker_of_each_bin.tar.gz",
        # unbinned=f"{outdir}/results/06_binning/maxbin2/{{sample_pool}}/{{sample_pool}}_bins/{{sample_pool}}_bin.noclass",
        # log=f"{outdir}/results/06_binning/maxbin2/{{sample_pool}}/{{sample_pool}}_bins/{{sample_pool}}_bin.log",
        # tooshort=f"{outdir}/results/06_binning/maxbin2/{{sample_pool}}/{{sample_pool}}_bins/{{sample_pool}}_bin.tooshort" #@Todo: make this temp.

    params:
        bin_prefix=f"{outdir}/results/06_binning/maxbin2/{{sample_pool}}/{{sample_pool}}_bins/{{sample_pool}}_bin"
    threads:
        config['threads']
    resources:
        mem_mb=config['max_mem']
    conda:
        "../envs/maxbin2.yaml"
    shell:
        """
	mkdir -p {output.bin_dir}
	run_MaxBin.pl -contig {input.assembly} -abund {input.coverage} -out {params.bin_prefix} -thread {threads}
	"""


def resolve_tool_bins(wildcards):
    """
    Returns the list of bin files for a given tool.
    """
    if wildcards.tool == "metabat2":
        return checkpoints.metabat2.get(sample_pool=wildcards.sample_pool).output.bin_dir
    elif wildcards.tool == "maxbin2":
        return checkpoints.maxbin2.get(sample_pool=wildcards.sample_pool).output.bin_dir
    else:
        raise ValueError(f"Unknown tool: {wildcards.tool}")


rule dastool_contigs2bin:
    input:
        bins_dir=resolve_tool_bins  # Dynamically resolve the correct bin directory,
    output:
        tsv=f"{outdir}/results/06_binning/{{tool}}/{{sample_pool}}/{{sample_pool}}_contigs2bin.tsv"
    params:
        script=workflow.source_path("../scripts/Fasta_to_Contig2Bin.sh"),
        extension=lambda wildcards: binning_tools[wildcards.tool]
    shell:
        "bash {params.script} -e {params.extension} -i {input.bins_dir} > {output.tsv} "


def format_dastool_input(input_files):
    """
    Converts a list of input Contigs2Bin.tsv files into a comma-separated string.
    """
    return ",".join(input_files)


checkpoint dastool:
    input:
#        assembly=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta",
        assembly=f"{outdir}/results/03_assembly/size_filtered/{{sample_pool}}_{minsize}/contigs_{{sample_pool}}_{minsize}.fasta",
        # contigs2bin_files=expand(f"{outdir}/results/06_binning/{{tool}}/{{sample_pool}}/{{sample_pool}}_contigs2bin.tsv",
        #                     tool=binning_tools.keys(),
        #                     sample_pool=samples["sample_pool"]),
        contigs2bin_metabat=f"{outdir}/results/06_binning/metabat2/{{sample_pool}}/{{sample_pool}}_contigs2bin.tsv",
        contigs2bin_maxbin=f"{outdir}/results/06_binning/maxbin2/{{sample_pool}}/{{sample_pool}}_contigs2bin.tsv"
    output:
        # dastool_output=directory(f"{outdir}/results/06_binning/dastool/{{sample_pool}}"),
        # quality_reports=f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool_summary.tsv",
        bin_dir=directory(f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool_bins"),
        # c2bin=f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool_contigs2bin.tsv"

    params:
        dastool_output=f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}",
        input_list=lambda wildcards, input: format_dastool_input([input.contigs2bin_metabat, input.contigs2bin_maxbin])
        # input_list=lambda wildcards, input: format_dastool_input(input.contigs2bin_files) #Comma-separated list
    threads:
        config['threads']
    conda:
        "../envs/dastool.yaml"
    shell:
        "DAS_Tool -i {params.input_list} -l metabat2,maxbin -c {input.assembly} -o {params.dastool_output} -t {threads} --write_bins"


rule checkm2:
    input:
        dastool_dir=f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool_bins"
    output:
        checkm_output=f"{outdir}/results/06_binning/checkm2/{{sample_pool}}/quality_report.tsv",
        diamond_output=f"{outdir}/results/06_binning/checkm2/{{sample_pool}}/diamond_output/DIAMOND_RESULTS.tsv",
        protein_files=directory(f"{outdir}/results/06_binning/checkm2/{{sample_pool}}/protein_files")
    params:
        output_dir=f"{outdir}/results/06_binning/checkm2/{{sample_pool}}",
        db_path=config['checkm_db']
    threads:
        config['threads']
    resources:
        mem_mb=config['max_mem']
    log: f"{outdir}/logs/{{sample_pool}}_checkm2.log"
    conda:
        "../envs/checkm2.yaml"
    shell:
        "checkm2 predict --threads {threads} -x fa --input {input.dastool_dir} --output-directory {params.output_dir} --force --database_path {params.db_path} 2> {log}"

rule BAT:
    input:
        # drep_dir=f"{outdir}/results/06_binning/drep/dereplicated_genomes",
        dastool_dir=f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool_bins",
        proteins= {rules.prodigal.output.faa},
        alignment= f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/{{sample_pool}}.alignment.diamond"
    output:
        bat_class=f"{outdir}/results/06_binning/BAT/{{sample_pool}}/{{sample_pool}}.bin2classification.txt",
        bat_names=f"{outdir}/results/06_binning/BAT/{{sample_pool}}/{{sample_pool}}.bin2classification.names.txt",
        bat_summary=f"{outdir}/results/06_binning/BAT/{{sample_pool}}/{{sample_pool}}.bin2classification.names.summarise.txt"
    params:
        output_dir=f"{outdir}/results/06_binning/BAT/{{sample_pool}}",
        db_path=config['CAT_database'],
        tax_path=config['CAT_taxonomy'],
        prefix=f"{{sample_pool}}"
    threads:
        config['threads']
    resources:
        mem_mb=config['max_mem']
    log: f"{outdir}/logs/{{sample_pool}}_bat.log"
    conda:
        "../envs/cat.yaml"
    shell:
        """ 
        mkdir -p {params.output_dir}
        CAT_pack bins -b {input.dastool_dir} -d {params.db_path} -t {params.tax_path} -p {input.proteins} \
         -a {input.alignment} -n {threads} -o {params.output_dir}{params.prefix} 2> {log}
        CAT_pack add_names -i {output.bat_class} -o {output.bat_names} -t {params.tax_path} --only_official --exclude_scores
        CAT_pack summarise -i {output.bat_names} -o {output.bat_summary}
        """

rule checkm2_to_drep_format:
    input:
        checkm2_output=f"{outdir}/results/06_binning/checkm2/{{sample_pool}}/quality_report.tsv"
    output:
        genome_info=f"{outdir}/results/06_binning/drep/checkm2_genomeinfo/{{sample_pool}}_genomeinfo.tsv"
    params:
        input_separator="\t",
        output_separator=",",
        genome_extension=".fa"
    run:
        # Read TSV file
        try:
            df = pd.read_csv(input.checkm2_output, sep=params.input_separator)
        except pd.errors.EmptyDataError:
            raise ValueError(f"Input file {input.checkm2_output} is empty")

        # Verify required columns exist
        required_cols = ['Name', 'Completeness', 'Contamination']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")

        # Add extension to genome names
        df['Name'] = df['Name'] + params.genome_extension

        # Rename columns for clarity
        column_mapping = {
            'Name': 'genome',
            'Completeness': 'completeness',
            'Contamination': 'contamination'
        }
        df = df.rename(columns=column_mapping)

        # Ensure only required columns are present and in correct order
        df = df[['genome', 'completeness', 'contamination']]

        # Save to CSV
        df.to_csv(output.genome_info, sep=params.output_separator, index=False)


rule combine_genome_info:
    input:
        genome_info_files=expand(f"{outdir}/results/06_binning/drep/checkm2_genomeinfo/{{sample_pool}}_genomeinfo.tsv",
            sample_pool=sorted(set(samples["sample_pool"])))
    output:
        combined_info=f"{outdir}/results/06_binning/drep/combined_genomeinfo.tsv"
    run:
        # Read and combine all genome info files
        dfs = []
        for file in input.genome_info_files:
            df = pd.read_csv(file,sep=",")  # Using comma as it's the output separator from previous rule
            dfs.append(df)

        # Concatenate all dataframes
        combined_df = pd.concat(dfs,ignore_index=True)

        # Save combined dataframe
        combined_df.to_csv(output.combined_info,sep=",",index=False)


checkpoint dereplicate_bins:
    input:
        bins_dir = expand(f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool_bins",
                  sample_pool=sorted(set(samples["sample_pool"]))
                  ),
        combined_info=f"{outdir}/results/06_binning/drep/combined_genomeinfo.tsv"
    output:
        dereplicated_bins = directory(f"{outdir}/results/06_binning/drep/dereplicated_genomes")
    params:
        drep_output = f"{outdir}/results/06_binning/drep",
        bin_dirs = lambda wildcards, input: ' '.join([f"{dir}/*.fa" for dir in sorted(set(input.bins_dir))])
    threads:
        config['threads']
    resources:
        mem_mb=config['max_mem']
    log:
        debug_log = f"{outdir}/results/06_binning/drep/drep_rule.log"
    conda:
        "../envs/drep.yaml"
    shell:
        """
        if [ $(echo "{input.bins_dir}" | tr ' ' '\n' | wc -l) -gt 1 ]; then
            echo "bin_dirs: {params.bin_dirs}" >> {log.debug_log}
            dRep dereplicate {params.drep_output} -g {params.bin_dirs} -p {threads} --genomeInfo {input.combined_info}
        else
            mkdir -p {output.dereplicated_bins}
            cp {input.bins_dir}/*.fa {output.dereplicated_bins}/
        fi
        """

