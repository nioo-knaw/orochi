# Script to calculate bin-level mean depths from contig-level coverage data. Based on the fairy-generated table.
# Each contig's depth is weighted by its length to compute accurate bin-level averages.
#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import re


def rename_and_concat_coverage_files(coverage_files, output_file):
    """
    For each coverage_<POOL>.tsv in coverage_files:
      - prepend <POOL>_ to contigName
      - concatenate all tables
      - write a single TSV to output_file

    Designed to be used as a single step in a Snakemake rule (no intermediate files required).
    """
    if not coverage_files:
        raise ValueError("coverage_files is empty")

    dfs = []
    expected_cols = None

    for coverage_file in coverage_files:
        base = os.path.basename(coverage_file)
        m = re.match(r"^coverage_(.+?)\.tsv$", base)
        if not m:
            raise ValueError(f"Expected filename like coverage_<POOL>.tsv, got: {base}")
        pool_id = m.group(1)

        df = pd.read_csv(coverage_file, sep="\t")
        if "contigName" not in df.columns:
            raise KeyError(f"Expected column 'contigName' in {coverage_file}, got columns: {list(df.columns)}")

        df["contigName"] = df["contigName"].astype(str).map(lambda x: f"{pool_id}_{x}")

        # Ensure all inputs have identical columns (helps catch mismatched sample columns early)
        cols = list(df.columns)
        if expected_cols is None:
            expected_cols = cols
        elif cols != expected_cols:
            raise ValueError(
                "Column mismatch while concatenating coverage files.\n"
                f"First file columns: {expected_cols}\n"
                f"{coverage_file} columns: {cols}"
            )

        dfs.append(df)

    out = pd.concat(dfs, axis=0, ignore_index=True)
    out.to_csv(output_file, sep="\t", index=False)


def rename_and_concat_contig2bin_files(contig2bin_files, output_file, write_header=False):
    """
    For each <POOL>_DASTool_contig2bin.tsv in contig2bin_files:
      - prepend <POOL>_ to contig IDs (col 1)
      - concatenate all files
      - write to output_file

    Input format (no header):
        contig<TAB>bin

    Example filename:
        A_DASTool_contig2bin.tsv  -> pool_id = "A"
    """
    if not contig2bin_files:
        raise ValueError("contig2bin_files is empty")

    dfs = []
    for path in contig2bin_files:
        base = os.path.basename(path)
        m = re.match(r"^(.+?)_DASTool_contig2bin\.tsv$", base)
        if not m:
            raise ValueError(f"Expected filename like <POOL>_DASTool_contig2bin.tsv, got: {base}")
        pool_id = m.group(1)

        df = pd.read_csv(path, sep="\t", header=None, names=["contigName", "bin"], dtype=str)
        if df.shape[1] != 2:
            raise ValueError(f"Expected exactly 2 columns in {path}, got {df.shape[1]}")

        df["contigName"] = df["contigName"].map(lambda x: f"{pool_id}_{x}")
        dfs.append(df)

    out = pd.concat(dfs, axis=0, ignore_index=True)

    out.to_csv(
        output_file,
        sep="\t",
        index=False,
        header=write_header,
    )


# Function to remove the contigs that are not in a bin from the coverage table
def filter_coverage_by_binned_contigs(coverage_file, mapping_file, output_file):
    """
    Filter a coverage TSV to keep only contigs that appear in the mapping file,
    preserving the contig order from the mapping file, and add the bin column
    (as the 2nd column) from the mapping.

    Assumptions:
      - coverage_file has a column 'contigName' plus coverage/depth columns
      - mapping_file has columns 'contigName' and 'bin' (or is a 2-col TSV without header)
      - mapping_file contains ONLY contigs that are assigned to bins

    Output:
      - columns: contigName, bin, <all original coverage columns except any existing 'bin'>
      - row order follows the mapping file's contig order
    """
    cov = pd.read_csv(coverage_file, sep="\t")
    if "contigName" not in cov.columns:
        raise KeyError(
            f"Expected column 'contigName' in coverage file {coverage_file}, got columns: {list(cov.columns)}"
        )

    # Read mapping with or without header; normalize to ['contigName','bin']
    map_df = pd.read_csv(mapping_file, sep="\t", dtype=str)
    if set(map_df.columns) >= {"contigName", "bin"}:
        map_df = map_df[["contigName", "bin"]]
    else:
        map_df = pd.read_csv(mapping_file, sep="\t", header=None, names=["contigName", "bin"], dtype=str)

    map_df["contigName"] = map_df["contigName"].astype(str)
    map_df["bin"] = map_df["bin"].astype(str)

    # Preserve mapping order; keep first bin if duplicates exist
    map_df = map_df.drop_duplicates(subset=["contigName"], keep="first")
    ordered_contigs = map_df["contigName"].tolist()

    # Index coverage by contigName for fast ordered lookup
    cov_idx = cov.copy()
    cov_idx["contigName"] = cov_idx["contigName"].astype(str)

    if "bin" in cov_idx.columns:
        cov_idx = cov_idx.drop(columns=["bin"])

    cov_idx = cov_idx.set_index("contigName", drop=False)

    present_contigs = [c for c in ordered_contigs if c in cov_idx.index]
    out = cov_idx.loc[present_contigs].reset_index(drop=True)

    # Add bin column in mapping order (aligned to 'present_contigs')
    bin_map = map_df.set_index("contigName")["bin"]
    out.insert(1, "bin", [bin_map.loc[c] for c in present_contigs])

    out.to_csv(output_file, sep="\t", index=False)


def _simplify_sample_colname(colname):
    """
    Convert a coverage table column name that may be a full path to a simple sample ID.
    Examples:
      /path/to/S1_filt_1.fastq.gz        -> S1
      /path/to/S1_filt_1.fastq.gz-var    -> S1-var
      S2_filt_1.fastq.gz                 -> S2
    If it doesn't look like a path/fastq column, returns the original name (collapsed to first '_' chunk).
    """
    s = str(colname)

    is_var = s.endswith("-var")
    if is_var:
        s = s[:-4]

    base = os.path.basename(s)

    # Usually looks like S1_filt_1.fastq.gz -> take "S1"
    sample = base.split("_")[0]
    return f"{sample}-var" if is_var else sample

def _suffix_sample_columns(df, suffix, exclude_cols=None):
    """
    Rename columns that look like sample IDs by appending a suffix.
    Example: S1 -> S1_meanDepth

    exclude_cols: columns that should never be suffixed (e.g. bin, binLen).
    """
    if exclude_cols is None:
        exclude_cols = set()
    else:
        exclude_cols = set(exclude_cols)

    rename_map = {}
    for c in df.columns:
        if c in exclude_cols:
            continue
        rename_map[c] = f"{c}_{suffix}"
    return df.rename(columns=rename_map)

def make_bin_coverage_from_binned_contig_coverage(binned_coverage_file, output_file, sample_value_suffix="meanDepth"):
    """
    Compute per-bin (length-weighted) mean depth from a *binned contig coverage* table.

    Output:
      - one row per bin
      - columns: bin, binLen, <sample columns like S1_meanDepth, S2_meanDepth, ...>
    """
    df = pd.read_csv(binned_coverage_file, sep="\t")

    required = {"contigName", "bin", "contigLen"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns in {binned_coverage_file}: {sorted(missing)}")

    meta_cols = {"contigName", "bin", "contigLen", "totalAvgDepth"}
    depth_cols = [c for c in df.columns if c not in meta_cols and not str(c).endswith("-var")]

    if not depth_cols:
        raise ValueError(
            "No depth columns detected. Expected depth columns besides contigName/bin/contigLen/totalAvgDepth "
            "and excluding *-var columns."
        )

    df["contigLen"] = pd.to_numeric(df["contigLen"], errors="raise")
    for c in depth_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    weighted = df[depth_cols].mul(df["contigLen"], axis=0)
    bin_weighted_sum = weighted.groupby(df["bin"]).sum()
    bin_len_sum = df.groupby("bin")["contigLen"].sum()
    bin_depth = bin_weighted_sum.div(bin_len_sum, axis=0)

    # Simplify sample column names (paths -> S1/S2/...)
    bin_depth = bin_depth.rename(columns=_simplify_sample_colname)

    # Add suffix to indicate what the numbers represent
    bin_depth = _suffix_sample_columns(bin_depth, sample_value_suffix)

    # Add bin length and write output
    bin_depth.insert(0, "binLen", bin_len_sum)
    bin_depth.index.name = "bin"
    bin_depth.reset_index().to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--coverage",
        nargs="+",
        required=True,
        help="One or more contig coverage tables (TSV).",
    )
    parser.add_argument(
        "-m",
        "--mapping",
        nargs="+",
        required=True,
        help="One or more contig-to-bin mapping files (TSV).",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory path.",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    tmp_cov = os.path.join(args.outdir, "temp_renamed_concat_coverage.tsv")
    tmp_map = os.path.join(args.outdir, "temp_renamed_concat_contig2bin.tsv")
    binned_cov = os.path.join(args.outdir, "binned_only_coverage.tsv")
    mag_depth_out = os.path.join(args.outdir, "mag_depth.tsv")

    rename_and_concat_coverage_files(coverage_files=args.coverage, output_file=tmp_cov)
    rename_and_concat_contig2bin_files(contig2bin_files=args.mapping, output_file=tmp_map, write_header=True)

    filter_coverage_by_binned_contigs(
        coverage_file=tmp_cov,
        mapping_file=tmp_map,
        output_file=binned_cov,
    )

    make_bin_coverage_from_binned_contig_coverage(
        binned_coverage_file=binned_cov,
        output_file=mag_depth_out,
    )

    os.remove(tmp_cov)
    os.remove(tmp_map)

    print(f"Bin-level mean depths calculated and saved to: {mag_depth_out}")
