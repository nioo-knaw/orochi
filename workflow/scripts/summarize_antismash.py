#!/usr/bin/env python3

"""Summarise antiSMASH GenBank output into a TSV table."""

from pathlib import Path
import csv

from Bio import SeqIO


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


def summarize_genbank_files(antismash_dir, taxon, output_tsv):
    """Create a TSV summary from antiSMASH GenBank region files."""
    antismash_dir = Path(antismash_dir)
    output_tsv = Path(output_tsv)

    genbank_files = sorted(
        list(antismash_dir.glob("*.region*.gbk"))
        + list(antismash_dir.glob("*001.gbk"))
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
                products = [
                    get_region_product(feature)
                    for feature in region_features
                ]
            else:
                products = [""]

            for region_number, product in enumerate(products, start=1):
                rows.append({
                    "sample_pool": snakemake.wildcards.sample_pool,
                    "taxon": taxon,
                    "bgc_id": f"{seq_record.id}_region{region_number}",
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
    antismash_dir=snakemake.input.antismash_dir,
    taxon=snakemake.params.taxon,
    output_tsv=snakemake.output.tsv,
)