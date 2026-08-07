#!/usr/bin/env python3

"""Combine the bacterial and fungal antiSMASH summaries into one BGC table.

One row per BGC region, annotated with the CAT contig-level taxonomy and --
for bacterial BGCs on a binned contig -- the bin it belongs to. Fungal BGCs
are never binned, so their bin stays N/A.
"""

from pathlib import Path

import pandas as pd

# Rank columns added by CAT_pack add_names --only_official. Ordered from
# highest to lowest rank so the last non-empty one is the most specific
# classification available for a contig.
CAT_RANK_COLUMNS = [
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

# CAT writes these instead of a name when a rank could not be resolved.
CAT_MISSING_VALUES = {"no support", "NA", "not classified", ""}

OUTPUT_COLUMNS = [
    "sample_pool",
    "taxon",
    "bgc_id",
    "contig_id",
    "bin_id",
    "taxonomy",
    "taxonomy_lowest",
    "bgc_product",
    "bgc_start",
    "bgc_end",
    "nr_genes",
]


def contig_key(series):
    """Normalise contig identifiers so the different tools' names match.

    summarize_antismash.py stores the GenBank *description*, which can carry
    trailing text (e.g. "k141_123 flag=1 multi=3.0"), while contig2bin and
    CAT use the bare contig name. Comparing on the first whitespace-delimited
    token makes both merges work either way.
    """
    return (
        series.astype(str)
        .str.strip()
        .str.split(n=1)
        .str[0]
    )


def read_summary(path):
    """Read one antiSMASH summary TSV, tolerating an empty/absent file."""
    try:
        summary_df = pd.read_csv(path, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return pd.DataFrame()

    return summary_df


def load_bgc_summaries(bacterial_path, fungal_paths):
    """Stack the bacterial and fungal summaries into one frame.

    fungal_paths is a (possibly empty) list: fungiSMASH is only run when the
    assembly unit has eukaryotic contigs with predicted genes.
    """
    frames = []

    bacterial_df = read_summary(bacterial_path)

    if not bacterial_df.empty:
        frames.append(bacterial_df)

    for fungal_path in fungal_paths:
        fungal_df = read_summary(fungal_path)

        if not fungal_df.empty:
            frames.append(fungal_df)

    if not frames:
        return pd.DataFrame(
            columns=[
                "sample_pool",
                "taxon",
                "bgc_id",
                "contig_id",
                "bgc_product",
                "bgc_start",
                "bgc_end",
                "nr_genes",
            ]
        )

    return pd.concat(frames, ignore_index=True)


def load_contig2bin(path):
    """Read the DASTool contig-to-bin mapping (headerless, two columns)."""
    try:
        c2b_df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["contig_id", "bin_id"],
        )
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return pd.DataFrame(columns=["contig_key", "bin_id"])

    c2b_df["contig_key"] = contig_key(c2b_df["contig_id"])

    return c2b_df[["contig_key", "bin_id"]].drop_duplicates("contig_key")


def load_cat_taxonomy(path):
    """Read CAT contig classifications into contig_key/lineage columns.

    The header line itself starts with "#" (e.g. "# contig\tclassification
    \t..."), so comment="#" must NOT be used here -- it would drop the header
    and promote the first data row to column names.
    """
    empty = pd.DataFrame(
        columns=["contig_key", "taxonomy", "taxonomy_lowest"]
    )

    try:
        cat_df = pd.read_csv(path, sep="\t", dtype=str)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return empty

    if cat_df.empty:
        return empty

    cat_df = cat_df.rename(columns={cat_df.columns[0]: "contig_id"})
    cat_df["contig_key"] = contig_key(cat_df["contig_id"])

    # CAT_pack's own "lineage" column holds taxids ("1;131567;2;..."), not
    # readable names, so build the lineage from the named rank columns.
    present_ranks = [
        column for column in CAT_RANK_COLUMNS if column in cat_df.columns
    ]

    if not present_ranks:
        return empty

    def build_lineage(row):
        names = []

        for value in row:
            if pd.isna(value):
                continue

            value = str(value).strip()

            if value in CAT_MISSING_VALUES:
                continue

            names.append(value)

        return ";".join(names)

    cat_df["taxonomy"] = cat_df[present_ranks].apply(build_lineage, axis=1)

    # Lowest available rank, kept alongside the full lineage -- the full
    # string stays in the TSV and in the cell tooltip, while the table itself
    # only has room for the most specific name.
    cat_df["taxonomy_lowest"] = cat_df["taxonomy"].apply(
        lambda lineage: lineage.split(";")[-1] if lineage else "N/A"
    )

    return (
        cat_df[["contig_key", "taxonomy", "taxonomy_lowest"]]
        .drop_duplicates("contig_key")
    )


def build_table(bgc_df, c2b_df, cat_df, sample_pool):
    """Annotate the BGC rows with bin and taxonomy, one row per BGC."""
    if bgc_df.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    bgc_df = bgc_df.copy()
    bgc_df["contig_key"] = contig_key(bgc_df["contig_id"])

    merged = bgc_df.merge(c2b_df, on="contig_key", how="left")
    merged = merged.merge(cat_df, on="contig_key", how="left")

    if "sample_pool" not in merged.columns:
        merged["sample_pool"] = sample_pool
    else:
        merged["sample_pool"] = merged["sample_pool"].fillna(sample_pool)

    for column in ("bin_id", "taxonomy", "taxonomy_lowest", "bgc_product"):
        if column not in merged.columns:
            merged[column] = "N/A"

        merged[column] = (
            merged[column]
            .fillna("N/A")
            .replace("", "N/A")
        )

    if "nr_genes" in merged.columns:
        merged["nr_genes"] = (
            pd.to_numeric(merged["nr_genes"], errors="coerce")
            .fillna(0)
            .astype(int)
        )
    else:
        merged["nr_genes"] = 0

    for column in OUTPUT_COLUMNS:
        if column not in merged.columns:
            merged[column] = "N/A"

    return merged[OUTPUT_COLUMNS].sort_values(
        ["taxon", "bgc_id"], kind="stable"
    )


def main():
    bacterial_path = Path(snakemake.input.bacterial_summary)
    fungal_paths = [Path(path) for path in snakemake.input.fungal_summary]

    contig2bin_path = Path(snakemake.input.contig2bin)
    cat_path = Path(snakemake.input.cat_taxonomy)

    output_path = Path(snakemake.output.combined_table)
    sample_pool = str(snakemake.params.sample_pool)

    bgc_df = load_bgc_summaries(bacterial_path, fungal_paths)
    c2b_df = load_contig2bin(contig2bin_path)
    cat_df = load_cat_taxonomy(cat_path)

    table = build_table(bgc_df, c2b_df, cat_df, sample_pool)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, sep="\t", index=False)

    n_bacterial = int((table["taxon"] == "bacteria").sum())
    n_fungal = int((table["taxon"] == "fungi").sum())

    print(f"Sample pool: {sample_pool}")
    print(f"Combined BGC table: {output_path}")
    print(f"Total BGCs: {len(table)}")
    print(f"Bacterial BGCs: {n_bacterial}")
    print(f"Fungal BGCs: {n_fungal}")


if __name__ == "__main__":
    main()