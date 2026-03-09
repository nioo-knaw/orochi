#!/usr/bin/env python

import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os


def size_filter(contig_file, outdir, size, sample_name):
    """
    Filter contigs based on their size.
    """
    print("Filtering contigs based on size...")

    os.makedirs(outdir, exist_ok=True)

    size = int(size) # Make sure size is an integer

    # Make sure input file exists
    if not os.path.isfile(contig_file):
        sys.exit(f"Error: Input file '{contig_file}' not found.")

    file_name = f"contigs_{sample_name}_{size}.fasta"
    total_contigs = 0
    kept_contigs = 0
    removed_contigs = 0

    with open(contig_file, 'r') as in_handle, open(os.path.join(outdir, file_name), 'w') as filtered_contigs:
        to_write = []

        for title, seq in SimpleFastaParser(in_handle):
            total_contigs += 1
            if len(seq) >= int(size):
                to_write.append(f">{title}\n{seq}\n")
                kept_contigs += 1
            else:
                removed_contigs += 1

        if total_contigs == 0:
            sys.exit(f"Error: No contigs found in '{contig_file}'. Check input.")

        print(f"Writing sequences > {size} to {filtered_contigs}")
        filtered_contigs.writelines(to_write) # Write all sequences > size at once


    print(f"Total contigs checked:\t{total_contigs}")
    print(f"Number of contigs >= {size} bp:\t{kept_contigs}")
    print(f"Number of contigs < {size} bp and removed:\t{removed_contigs}")

    if kept_contigs == 0:
        sys.exit(f"No contigs longer than {size} bp detected.\n")
    if removed_contigs == 0:
        pass  # This can return a flag for skipping downstream steps.


if __name__ == "__main__":
    # Parse arguments from Snakemake
    contig_file = sys.argv[2]
    outdir = sys.argv[4]
    size = int(sys.argv[6])
    sample_name = sys.argv[8]

    size_filter(contig_file, outdir, size, sample_name)
