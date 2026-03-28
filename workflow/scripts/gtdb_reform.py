#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse


def process_abundance_table(input_file, output_file):
    """
    Process the merged abundance table:
    1. Replace all ';' with '|'
    2. Find all rows containing species level (with 's__')
    3. Copy those rows to the end of the file and append '|t__' to the clade_name
    4. Save the result to the output file
    """
    # Read all lines from input file
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        return
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    # Step 1: Replace all ';' with '|'
    processed_lines = [line.replace(';', '|') for line in lines]

    # Step 2: Collect all species-level rows (those containing 's__' in clade_name)
    species_rows = []
    for line in processed_lines:
        line = line.rstrip('\n')
        if not line.strip() or line.startswith('#'):
            continue
        if '\t' not in line:
            continue
        parts = line.split('\t')
        clade_name = parts[0]
        # Check if this is a species-level entry
        if 's__' in clade_name:
            species_rows.append(line + '\n')

    # Step 3: Create new rows with '|t__' appended to clade_name
    new_species_rows = []
    for row in species_rows:
        parts = row.rstrip('\n').split('\t')
        original_clade = parts[0]
        # Append |t__ to the first column
        new_clade = original_clade + '|t__'
        # Rebuild the row
        new_row = new_clade + '\t' + '\t'.join(parts[1:]) + '\n'
        new_species_rows.append(new_row)

    # Step 4: Combine original processed lines + new species rows
    final_lines = processed_lines + new_species_rows

    # Step 5: Write to output file
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.writelines(final_lines)
        print(f"Processing completed.")
        # print(f"Output saved to: {output_file}")
        # print(f"Added {len(new_species_rows)} new rows with |t__")
    except Exception as e:
        print(f"Error writing output file: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Reformat merged abundance table: replace ; with |, "
                    "and append |t__ rows for species-level entries."
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Path to input file (merged_abundance_table.txt)')
    parser.add_argument('-o', '--output', required=True,
                        help='Path to output file (e.g. processed_abundance_table.txt)')

    args = parser.parse_args()

    process_abundance_table(args.input, args.output)


if __name__ == "__main__":
    main()