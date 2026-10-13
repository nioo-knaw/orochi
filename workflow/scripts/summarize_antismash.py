#!/usr/bin/env python3

"""Summarise antiSMASH GenBank output into a TSV table."""

from pathlib import Path
import csv
import re
import sys
from Bio import SeqIO

log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log


def get_structured_antismash_value(seq_record, key, default=""):
    """Safely fetch values from the antiSMASH structured comment block."""
    return (
        seq_record.annotations
        .get("structured_comment", {})
        .get("antiSMASH-Data", {})
        .get(key, default)
    )


def get_region_product(feature):
    """Return antiSMASH region product as a comma-separated string."""
    products = feature.qualifiers.get("product", [])
    if isinstance(products, list):
        return ",".join(products)
    return str(products)


def get_region_number(feature, gbk_file, fallback):
    """Return antiSMASH's own number for a region.

    Each *.regionNNN.gbk holds a single region, so the region's position
    within the file is always 1 and cannot be used to tell the regions of one
    contig apart. antiSMASH's own numbering lives in the region feature's
    /region_number qualifier, with the file name (".region002.gbk") as a
    fallback for output that lacks the qualifier.
    """
    numbers = feature.qualifiers.get("region_number", [])

    if numbers:
        try:
            return int(str(numbers[0]).strip())
        except ValueError:
            pass

    match = re.search(r"\.region(\d+)\.gbk$", gbk_file.name)

    if match:
        return int(match.group(1))

    return fallback


def summarize_genbank_files(antismash_dir, taxon, output_tsv):
    """Create a TSV summary from antiSMASH GenBank region files."""
    antismash_dir = Path(antismash_dir)
    output_tsv = Path(output_tsv)

    genbank_files = sorted(
        list(antismash_dir.glob("*.region*.gbk"))
    )

    rows = []

    for gbk_file in genbank_files:
        for seq_record in SeqIO.parse(gbk_file, "genbank"):
            cds_count = sum(
                1 for feature in seq_record.features
                if feature.type == "CDS"
            )

            region_features = [
                feature for feature in seq_record.features
                if feature.type == "region"
            ]

            if region_features:
                regions = [
                    (
                        get_region_number(feature, gbk_file, position),
                        get_region_product(feature),
                    )
                    for position, feature in enumerate(region_features, start=1)
                ]
            else:
                regions = [(1, "")]

            for region_number, product in regions:
                rows.append({
                    "sample_pool": snakemake.wildcards.sample_pool,
                    "taxon": taxon,
                    "bgc_id": f"{seq_record.id}_region{region_number}",
                    "region_number": region_number,
                    "contig_id": seq_record.description,
                    "bgc_product": product,
                    "bgc_start": get_structured_antismash_value(seq_record, "Orig. start"),
                    "bgc_end": get_structured_antismash_value(seq_record, "Orig. end"),
                    "nr_genes": cds_count,
                    "genbank_file": gbk_file.name,
                })

    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "sample_pool",
        "taxon",
        "bgc_id",
        "region_number",
        "contig_id",
        "bgc_product",
        "bgc_start",
        "bgc_end",
        "nr_genes",
        "genbank_file",
    ]

    with output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


summarize_genbank_files(
    antismash_dir=snakemake.params.antismash_dir,
    taxon=snakemake.params.taxon,
    output_tsv=snakemake.output.tsv,
)