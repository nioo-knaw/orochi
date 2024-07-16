def get_extension():
    if config['normalize_reads'] == 'Yes':
        sample_extension = "_normalized_"
    elif config['normalize_reads'] == "No":
        sample_extension = "_filt_"
    return sample_extension

def get_forward_files(wildcards, extension):
    return ["results/02_filtered_reads/" + sample + extension + "1.fastq.gz" for sample in
            samples[samples["sample_pool"] == wildcards.sample_pool]["sample"].values]


def get_rev_files(wildcards, extension):
    return ["results/02_filtered_reads/" + sample + extension + "2.fastq.gz" for sample in
            samples[samples["sample_pool"] == wildcards.sample_pool]["sample"].values]


# def assembly_input():
#     # Determine the normalization status based on the config
#     norm_status = config['normalize_reads']
#     print(f"Normalization: {norm_status}")
#
#     # Construct the file paths for both forward and reverse reads
#     if norm_status == "Yes":
#         return [
#             "results/03_assembly/coassembly/pools/{sample_pool}_normalized_f.fastq",  # Forward read
#             "results/03_assembly/coassembly/pools/{sample_pool}_normalized_r.fastq" # Reverse read
#         ]
#     else:
#         return [
#             "results/03_assembly/coassembly/pools/{sample_pool}_forward.fastq",  # Forward read
#             "results/03_assembly/coassembly/pools/{sample_pool}_rev.fastq"  # Reverse read
#         ]
