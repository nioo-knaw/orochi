#!/usr/bin/env python3
"""
Filter antiSMASH GenBank output files to only include contigs from a specific bin.
"""

import sys
from pathlib import Path
from Bio import SeqIO


def filter_genbank_by_bin(antismash_dir, contig2bin_file, bin_id, output_dir):
    """
    Extract GenBank records for contigs belonging to a specific bin.
    
    Args:
        antismash_dir: Path to antiSMASH output directory with GenBank files
        contig2bin_file: Path to contig2bin TSV file
        bin_id: Bin ID to filter for
        output_dir: Output directory for filtered GenBank files
    """
    antismash_dir = Path(antismash_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load contig to bin mapping
    contig_to_bin = {}
    with open(contig2bin_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    contig_id, bin_name = parts[0], parts[1]
                    contig_to_bin[contig_id] = bin_name
    
    # Get contigs for this bin
    bin_contigs = {contig for contig, bin_name in contig_to_bin.items() 
                   if bin_name == bin_id}
    
    if not bin_contigs:
        print(f"Warning: No contigs found for bin {bin_id}", file=sys.stderr)
        return 0
    
    # Find all GenBank files
    genbank_files = sorted(
        list(antismash_dir.glob("*.region*.gbk")) + 
        list(antismash_dir.glob("*[0-9].gbk"))
    )
    
    if not genbank_files:
        print(f"Warning: No GenBank files found in {antismash_dir}", file=sys.stderr)
        return 0
    
    regions_written = 0
    
    # Filter GenBank files
    for gbk_file in genbank_files:
        records_to_write = []
        
        try:
            for seq_record in SeqIO.parse(gbk_file, "genbank"):
                # Check if this contig belongs to the bin
                # antiSMASH sometimes modifies contig IDs, so we check if any bin contig
                # is a substring or if the record ID/description matches
                contig_match = False
                
                for bin_contig in bin_contigs:
                    if (bin_contig in seq_record.id or 
                        bin_contig in seq_record.description or
                        seq_record.id in bin_contig):
                        contig_match = True
                        break
                
                if contig_match:
                    records_to_write.append(seq_record)
        except Exception as e:
            print(f"Error parsing {gbk_file}: {e}", file=sys.stderr)
            continue
        
        # Write filtered records
        if records_to_write:
            output_file = output_dir / gbk_file.name
            SeqIO.write(records_to_write, output_file, "genbank")
            regions_written += len(records_to_write)
    
    return regions_written


# When called as a Snakemake script, parameters are available via snakemake object
if __name__ == "__main__":
    # Check if running as a Snakemake script
    try:
        # Access Snakemake variables
        antismash_dir = snakemake.params.antismash_dir
        contig2bin = snakemake.input.contig2bin
        bin_id = snakemake.params.bin_id
        output_dir = snakemake.output.filtered_dir
        
        n_regions = filter_genbank_by_bin(
            antismash_dir,
            contig2bin,
            bin_id,
            output_dir
        )
        
        print(f"Filtered {n_regions} BGC regions for bin {bin_id}")
        
    except NameError:
        # Not running as Snakemake script, use argparse
        import argparse
        
        parser = argparse.ArgumentParser(
            description="Filter antiSMASH GenBank files by bin"
        )
        parser.add_argument(
            "--antismash-dir",
            required=True,
            help="Path to antiSMASH output directory"
        )
        parser.add_argument(
            "--contig2bin",
            required=True,
            help="Path to contig2bin TSV file"
        )
        parser.add_argument(
            "--bin-id",
            required=True,
            help="Bin ID to filter for"
        )
        parser.add_argument(
            "--output-dir",
            required=True,
            help="Output directory for filtered GenBank files"
        )
        
        args = parser.parse_args()
        
        n_regions = filter_genbank_by_bin(
            args.antismash_dir,
            args.contig2bin,
            args.bin_id,
            args.output_dir
        )
        
        print(f"Filtered {n_regions} BGC regions for bin {args.bin_id}")