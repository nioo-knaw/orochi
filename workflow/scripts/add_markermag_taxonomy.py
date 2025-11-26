# Linking the phyloflash extracted SSU classifications to the markermag linkages by genome file
import pandas as pd
import os
from pathlib import Path
import csv
import argparse

def link_markermag_taxonomy(markermag_file, phyloflash_file, output_file):
    # Read the markermag file
    markermag_df = pd.read_csv(markermag_file, sep='\t')

    # Read the phyloflash file into a dictionary
    phyloflash_dict = {}
    with open(phyloflash_file) as file:
        tsv_file = csv.reader(file, delimiter=",")
        for line in tsv_file:
            if line[0].startswith("OTU"):
                continue
            otu_id = line[0]
            taxonomy = line[4]
            phyloflash_dict[otu_id] = taxonomy
            print(taxonomy)
        print(phyloflash_dict)

    # Prepare a list to hold the new data with taxonomy
    new_data = []

    # Iterate through the markermag dataframe and add taxonomy
    for index, row in markermag_df.iterrows():
        marker_id_raw = row['MarkerGene'] # Adjust column name as per your file
        marker_id = "_".join(marker_id_raw.split("_")[:2])
        taxonomy = phyloflash_dict.get(marker_id, 'unknown')
        new_row = row.tolist() + [taxonomy]
        new_data.append(new_row)
        print(new_row)

    # Write the new data to the output file
    with open(output_file, 'w', newline='') as tsvfile:
        taxwriter = csv.writer(tsvfile, delimiter='\t')
        # Write header
        header = list(markermag_df.columns) + ['taxonomy']
        taxwriter.writerow(header)
        # Write data rows
        for new_row in new_data:
            taxwriter.writerow(new_row)
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Link markermag taxonomy with phyloflash classifications')
    parser.add_argument('-m', '--markermag', help='Path to the markermag input file')
    parser.add_argument('-p', '--phyloflash', help='Path to the phyloflash input file')
    parser.add_argument('-o', '--output', help='Path to the output file')
    
    args = parser.parse_args()
    
    link_markermag_taxonomy(args.m, args.p, args.o)