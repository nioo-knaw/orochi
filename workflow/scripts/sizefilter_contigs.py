#!/usr/bin/env python

import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os


def size_filter(contig_file, outdir, size, sample_name):
    """
    Filter contigs based on their size.
    """
    print("Filtering contigs based on size...")
    file_name = f"contigs_{size}_{sample_name}.fasta"
    with open(contig_file) as in_handle:
        with open(os.path.join(outdir, file_name), 'w', newline='') as filtered_contigs:
            total_contigs = 0
            kept_contigs = 0
            removed_contigs = 0

            for title, seq in SimpleFastaParser(in_handle):
                total_contigs += 1
                if len(seq) >= int(size):
                    filtered_contigs.write(f">{title}\n{seq}\n")
                    kept_contigs += 1
                else:
                    removed_contigs += 1

            print(f"Total contigs checked:\t{total_contigs}")
            print(f"Number of contigs >= {size} bp:\t{kept_contigs}")
            print(f"Number of contigs < {size} bp and removed:\t{removed_contigs}")

            if kept_contigs == 0:
                sys.exit(f"No contigs longer than {size} bp detected.\n")
            if removed_contigs == 0:
                pass  # This can return a flag for skipping downstream steps.


if __name__ == "__main__":
    # Parse arguments from Snakemake
    contig_file = sys.argv[1]
    outdir = sys.argv[2]
    size = int(sys.argv[3])
    sample_name = sys.argv[4]

    size_filter(contig_file, outdir, size, sample_name)
