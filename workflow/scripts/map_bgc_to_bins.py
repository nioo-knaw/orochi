import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Map BGCs to MAGs and their taxonomy."
    )

    parser.add_argument("--bgc", required=True)
    parser.add_argument("--contig2bin", required=True)
    parser.add_argument("--bat", required=True)
    parser.add_argument("--markermag", required=True)
    parser.add_argument("--output", required=True)

    return parser.parse_args()

def main():
    args = parse_args()

    # Load BGC summary
    bgc_df = pd.read_csv(args.bgc, sep="\t")

    # Load contig to bin mapping
    c2b_df = pd.read_csv(
        args.contig2bin,
        sep="\t",
        header=None,
        names=["contig_id", "bin_id"]
    )

    # Load BAT taxonomy. The header line itself starts with "#" (e.g.
    # "# bin\tclassification\treason\t..."), so comment="#" must NOT be
    # used here -- it would drop the header and treat the first data row
    # as column names instead.
    bat_df = pd.read_csv(args.bat, sep="\t")
    bat_df = bat_df.rename(columns={bat_df.columns[0]: "bin_id"})
    # CAT_pack bins was run with -s .fa, but add_names does not strip
    # that suffix from the bin name (e.g. "A_bin.001.fa"), while
    # DASTool_contig2bin.tsv uses the bare bin id (e.g. "A_bin.001").
    # Strip it so the merge below on bin_id actually matches.
    bat_df["bin_id"] = bat_df["bin_id"].str.removesuffix(".fa")

    # CAT_pack's own "lineage" column is a taxid string (e.g.
    # "1;131567;2;..."), not a human-readable name. Build a readable
    # lineage from the named rank columns added by --only_official,
    # skipping ranks with no support/NA.
    rank_cols = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    bat_df["lineage"] = bat_df[rank_cols].apply(
        lambda row: ";".join(
            str(v) for v in row if pd.notna(v) and v not in ("no support", "NA")
        ),
        axis=1
    )
    # Lowest available rank, kept alongside the full lineage -- the full
    # string stays in the tsv for downstream figure scripts, while the
    # HTML report only has room to show the lowest rank per bin.
    bat_df["lineage_lowest"] = bat_df["lineage"].apply(
        lambda s: s.split(";")[-1] if s else "N/A"
    )

    # Merge BGC with bin assignment
    bgc_bins = bgc_df.merge(c2b_df,on="contig_id",how="left")

    # Merge with taxonomy
    bgc_bins_tax = bgc_bins.merge(
        bat_df[["bin_id", "lineage", "lineage_lowest"]],
        on="bin_id",
        how="left"
    )

    # Load markerMAG taxonomy (16S-based genome-marker linkage). Not every
    # bin is reached by a 16S linkage, so this is a left join and missing
    # values are expected -- shown alongside, not instead of, BAT taxonomy.
    try:
        markermag_df = pd.read_csv(args.markermag, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        markermag_df = pd.DataFrame(columns=["GenomicSeq", "taxonomy"])

    if not markermag_df.empty and "GenomicSeq" in markermag_df.columns:
        markermag_df = markermag_df.rename(
            columns={"GenomicSeq": "bin_id", "taxonomy": "markermag_taxonomy"}
        )
        # "unknown" is add_markermag_taxonomy.py's placeholder for marker
        # genes with no phyloFlash classification -- not informative here.
        markermag_df = markermag_df[
            markermag_df["markermag_taxonomy"].notna() &
            (markermag_df["markermag_taxonomy"] != "unknown") &
            (markermag_df["markermag_taxonomy"].str.strip() != "")
        ].copy()
        markermag_df["markermag_taxonomy_lowest"] = markermag_df["markermag_taxonomy"].apply(
            lambda s: s.split(";")[-1] if isinstance(s, str) and s else "N/A"
        )
        # A bin can have multiple marker gene linkages (e.g. several 16S
        # copies); collapse to one row per bin. "||" separates distinct
        # full lineages (each internally ";"-delimited by rank) so it
        # stays unambiguous; the lowest-rank names use ", " since they're
        # single tokens.
        markermag_tax = markermag_df.groupby("bin_id").agg(
            markermag_taxonomy=("markermag_taxonomy", lambda x: " || ".join(sorted(set(x)))),
            markermag_taxonomy_lowest=("markermag_taxonomy_lowest", lambda x: ", ".join(sorted(set(x))))
        ).reset_index()
    else:
        markermag_tax = pd.DataFrame(columns=["bin_id", "markermag_taxonomy", "markermag_taxonomy_lowest"])

    bgc_bins_tax = bgc_bins_tax.merge(markermag_tax, on="bin_id", how="left")

    print(f"BGC records: {len(bgc_df)}")
    print(f"Contig-to-bin records: {len(c2b_df)}")
    print(f"BAT records: {len(bat_df)}")
    print(f"markerMAG records: {len(markermag_df)}")

    # Add indicator for unbinned contigs
    bgc_bins_tax["bin_status"] = bgc_bins_tax["bin_id"].apply(
        lambda x: "binned" if pd.notna(x) else "unbinned"
    )

    # Save combined table
    bgc_bins_tax.to_csv(args.output,sep="\t",index=False)

if __name__ == "__main__":
    main()