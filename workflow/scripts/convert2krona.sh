#!/bin/bash

# Usage: ./convert_to_krona.sh input.tsv > krona_input.txt

# Skip header and process lines
tail -n +2 "$1" | while IFS=$'\t' read -r contig classification reason lineage scores sk phylum class order family genus species; do
  # Skip unclassified contigs
  if [[ "$classification" != "taxid assigned" ]]; then
    continue
  fi

  # Build taxonomy path (skip 'no support' entries)
  taxonomy=()
  for level in "$sk" "$phylum" "$class" "$order" "$family" "$genus" "$species"; do
    if [[ "$level" != "no support" && "$level" != "NA" && -n "$level" ]]; then
      taxonomy+=("$level")
    fi
  done

  # Output in Krona format: count \t taxonomy path
  if [[ ${#taxonomy[@]} -gt 0 ]]; then
    echo -e "1\t${taxonomy[*]}" | sed 's/ /;/g'
  fi
done
