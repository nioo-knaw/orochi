"""Common functions and rules used in the workflow"""

def get_forward_files(wildcards, extension):
    return [f"{outdir}/results/02_filtered_reads/" + sample + extension + "1.fastq.gz" for sample in
            samples[samples["sample_pool"] == wildcards.sample_pool]["sample"].values]


def get_rev_files(wildcards, extension):
    return [f"{outdir}/results/02_filtered_reads/" + sample + extension + "2.fastq.gz" for sample in
            samples[samples["sample_pool"] == wildcards.sample_pool]["sample"].values]

def prefix(assembly_method, samples):
    if assembly_method == "coassembly":
        return samples["sample_pool"].unique()
    elif assembly_type == "single_assembly":
        return samples["sample"].unique()

rule downstream_test:
    input:
        f"{outdir}/results/03_assembly/coassembly/assembly_{{sample}}/{{sample}}_assembly.done"

    output:
        test_file=f"{outdir}/results/05_test/{{sample}}/{{sample}}_test.txt"
    run:
        shell("touch {output.test_file}")
