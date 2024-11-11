"""Common functions and rules used in the workflow"""

def get_forward_files(wildcards, extension):
    return [f"{outdir}/results/02_filtered_reads/" + sample + extension + "1.fastq.gz" for sample in
            samples[samples["sample_pool"] == wildcards.sample_pool]["sample"].values]


def get_rev_files(wildcards, extension):
    return [f"{outdir}/results/02_filtered_reads/" + sample + extension + "2.fastq.gz" for sample in
            samples[samples["sample_pool"] == wildcards.sample_pool]["sample"].values]


def get_forward_files2(wildcards, extension):
    # Get list of samples for the given pool
    sample_list = [sample for sample, pool in sample_to_pool.items() if pool == wildcards.sample_pool]

    # Generate forward files only for the correct pool
    forward_files = [f"{outdir}/results/02_filtered_reads/{sample}{extension}1.fastq.gz" for sample in sample_list]

    return forward_files


def get_rev_files2(wildcards, extension):
    # Get list of samples for the given pool
    sample_list = [sample for sample, pool in sample_to_pool.items() if pool == wildcards.sample_pool]

    # Generate reverse files only for the correct pool
    rev_files = [f"{outdir}/results/02_filtered_reads/{sample}{extension}2.fastq.gz" for sample in sample_list]

    return rev_files


def prefix(assembly_method, samples):
    if assembly_method == "coassembly":
        return samples["sample_pool"].unique()
    elif assembly_type == "single_assembly":
        return samples["sample"].unique()

rule downstream_test:
    input:
        # f"{outdir}/results/06_binning/mapping/{{sample}}/{{sample}}_assembly.fasta.bwt",
        f"{outdir}/results/06_binning/metabat2/{{sample}}/{{sample}}_bins.tsv"

    output:
        test_file=f"{outdir}/results/05_test/{{sample}}/{{sample}}_test.txt"
    run:
        shell("echo {input}")
        shell("touch {output.test_file}")
