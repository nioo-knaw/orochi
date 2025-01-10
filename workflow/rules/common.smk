"""Common functions and rules used in the workflow"""
import os
from glob import glob

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

def get_metabat_bins(wildcards):
    # Retrieve checkpoint output
    checkpoint_output = checkpoints.metabat2.get(sample_pool=wildcards.sample_pool).output.bin_dir
    # Return the list of bin files
    return sorted(glob(f"{checkpoint_output}/{wildcards.sample_pool}_bin.*.fa"))

def get_maxbin_bins(wildcards):
    # Retrieve checkpoint output
    checkpoint_output = checkpoints.maxbin2.get(sample_pool=wildcards.sample_pool).output.bin_dir
    # Return the list of bin files
    return sorted(glob(f"{checkpoint_output}/{wildcards.sample_pool}_bin.*.fasta"))

def get_dastool_bins(wildcards):
    # Retrieve checkpoint output
    checkpoint_output = checkpoints.dastool.get(sample_pool=wildcards.sample_pool).output.bin_dir
    # Return the list of bin files
    return sorted(glob(f"{checkpoint_output}/{wildcards.sample_pool}_bin.*.fa"))

def get_drep_bins(wildcards):
    # Retrieve checkpoint output
    checkpoint_output = checkpoints.dereplicate_bins.get().output.dereplicated_bins
    # Return the list of bin files
    return sorted(glob(f"{checkpoint_output}/*.fa"))

