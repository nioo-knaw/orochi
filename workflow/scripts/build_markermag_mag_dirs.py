#!/usr/bin/env python3

import argparse
import os
import shutil
from pathlib import Path

import pandas as pd


FASTA_SUFFIXES = {".fa", ".fna", ".fasta"}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build per-sample-pool MAG directories for MarkerMAG from dRep results. "
            "For each original MAG starting with a sample_pool prefix, the script finds "
            "the final dRep representative MAG from the same secondary cluster."
        )
    )
    parser.add_argument(
        "--drep-dir",
        required=True,
        type=Path,
        help="dRep output directory, e.g. results/06_binning/drep"
    )
    parser.add_argument(
        "--input-bins",
        required=True,
        type=Path,
        help="Text file containing original MAG paths, one path per line"
    )
    parser.add_argument(
        "--sample-pool",
        required=True,
        help="Current sample_pool or assembly unit, e.g. S4"
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        type=Path,
        help="Output MAG directory for this sample_pool"
    )
    parser.add_argument(
        "--copy-mode",
        choices=["symlink", "copy"],
        default="symlink",
        help="Use symlink or copy for representative MAG files"
    )
    parser.add_argument(
        "--done",
        required=True,
        type=Path,
        help="Done file to create after successful completion"
    )
    return parser.parse_args()


def normalized_keys(path):
    """
    Return possible keys for matching dRep table entries and file names.
    """
    p = Path(str(path))
    return {
        str(path),
        p.name,
        p.stem,
    }


def starts_with_sample_pool(filename, sample_pool):
    """
    Match MAGs belonging to sample_pool.

    Prefer strict prefix sample_pool + '_' to avoid S1 matching S10.
    Fallback to sample_pool prefix for pipelines that do not use underscores.
    """
    name = Path(filename).name

    if name.startswith(sample_pool + "_"):
        return True

    if name.startswith(sample_pool + "."):
        return True

    if name.startswith(sample_pool + "-"):
        return True

    return False


def find_cluster_column(cdb):
    preferred = [
        "secondary_cluster",
        "secondary_cluster_id",
        "Secondary_cluster",
        "cluster_secondary",
    ]

    for col in preferred:
        if col in cdb.columns:
            return col

    candidates = [
        col for col in cdb.columns
        if "secondary" in col.lower() and "cluster" in col.lower()
    ]

    if candidates:
        return candidates[0]

    candidates = [
        col for col in cdb.columns
        if "cluster" in col.lower()
    ]

    if candidates:
        return candidates[0]

    raise ValueError(
        "Could not find a secondary cluster column in Cdb.csv. "
        f"Available columns: {list(cdb.columns)}"
    )


def build_cdb_lookup(cdb):
    """
    Build lookup from genome path/name/stem to Cdb row index.
    """
    if "genome" not in cdb.columns:
        raise ValueError(
            "Cdb.csv does not contain a 'genome' column. "
            f"Available columns: {list(cdb.columns)}"
        )

    lookup = {}

    for idx, genome in cdb["genome"].items():
        for key in normalized_keys(genome):
            lookup[key] = idx

    return lookup


def find_cdb_index_for_path(path, cdb_lookup):
    for key in normalized_keys(path):
        if key in cdb_lookup:
            return cdb_lookup[key]
    return None


def collect_representatives(derep_genomes_dir):
    reps = []

    for file in sorted(derep_genomes_dir.iterdir()):
        if file.is_file() and file.suffix in FASTA_SUFFIXES:
            reps.append(file)

    return reps


def link_or_copy(src, dst, copy_mode):
    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists() or dst.is_symlink():
        dst.unlink()

    if copy_mode == "symlink":
        os.symlink(src.resolve(), dst)
    else:
        shutil.copy2(src, dst)


def main():
    args = parse_args()

    drep_dir = args.drep_dir
    derep_genomes_dir = drep_dir / "dereplicated_genomes"
    cdb_file = drep_dir / "data_tables" / "Cdb.csv"

    if not derep_genomes_dir.exists():
        raise FileNotFoundError(f"dRep dereplicated_genomes directory not found: {derep_genomes_dir}")

    if not args.input_bins.exists():
        raise FileNotFoundError(f"Input bins file not found: {args.input_bins}")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    original_bins = [
        Path(line.strip())
        for line in args.input_bins.read_text().splitlines()
        if line.strip()
    ]

    sample_bins = [
        path for path in original_bins
        if starts_with_sample_pool(path.name, args.sample_pool)
    ]

    representatives = collect_representatives(derep_genomes_dir)

    if not representatives:
        raise FileNotFoundError(f"No representative MAG files found in {derep_genomes_dir}")

    rep_by_key = {}
    for rep in representatives:
        for key in normalized_keys(rep):
            rep_by_key[key] = rep

    selected_reps = {}

    # Case 1: normal dRep output with Cdb.csv
    if cdb_file.exists():
        cdb = pd.read_csv(cdb_file)
        cluster_col = find_cluster_column(cdb)
        cdb_lookup = build_cdb_lookup(cdb)

        # Map cluster -> representative MAGs found in dereplicated_genomes
        cluster_to_reps = {}

        for rep in representatives:
            idx = find_cdb_index_for_path(rep, cdb_lookup)
            if idx is None:
                # Some dRep versions store only original input basename/path.
                # Try matching representative by name/stem through Cdb genome column.
                continue

            cluster = cdb.loc[idx, cluster_col]
            cluster_to_reps.setdefault(cluster, []).append(rep)

        for original_bin in sample_bins:
            idx = find_cdb_index_for_path(original_bin, cdb_lookup)

            if idx is None:
                # Fallback: if this original MAG itself survived dRep, use it.
                rep = None
                for key in normalized_keys(original_bin):
                    if key in rep_by_key:
                        rep = rep_by_key[key]
                        break

                if rep is not None:
                    selected_reps[rep.name] = rep
                else:
                    print(
                        f"WARNING: {original_bin} was not found in Cdb.csv and no representative was found. Skipped."
                    )
                continue

            cluster = cdb.loc[idx, cluster_col]
            reps_in_cluster = cluster_to_reps.get(cluster, [])

            if reps_in_cluster:
                # Usually one winner per secondary cluster.
                rep = sorted(reps_in_cluster, key=lambda x: x.name)[0]
                selected_reps[rep.name] = rep
            else:
                # Fallback: if the original MAG itself survived dRep.
                rep = None
                for key in normalized_keys(original_bin):
                    if key in rep_by_key:
                        rep = rep_by_key[key]
                        break

                if rep is not None:
                    selected_reps[rep.name] = rep
                else:
                    print(
                        f"WARNING: No representative found for {original_bin} in cluster {cluster}. Skipped."
                    )

    # Case 2: no Cdb.csv, for example when only one bin was copied manually
    else:
        print(f"WARNING: {cdb_file} not found. Falling back to prefix-based representative selection.")

        for rep in representatives:
            if starts_with_sample_pool(rep.name, args.sample_pool):
                selected_reps[rep.name] = rep

    # Clean old files in output directory
    for old in args.out_dir.iterdir():
        if old.is_file() or old.is_symlink():
            old.unlink()

    for rep_name, rep_path in sorted(selected_reps.items()):
        dst = args.out_dir / rep_name
        link_or_copy(rep_path, dst, args.copy_mode)

    if not selected_reps:
        print(
            f"WARNING: No representative MAGs selected for sample_pool={args.sample_pool}. "
            "MarkerMAG may fail if this sample has no MAGs."
        )

    args.done.parent.mkdir(parents=True, exist_ok=True)
    args.done.write_text(
        f"sample_pool\t{args.sample_pool}\n"
        f"n_original_bins\t{len(sample_bins)}\n"
        f"n_representative_mags\t{len(selected_reps)}\n"
    )

    print(f"Built MAG directory for {args.sample_pool}: {args.out_dir}")
    print(f"Original MAGs from this sample_pool: {len(sample_bins)}")
    print(f"Representative MAGs selected: {len(selected_reps)}")


if __name__ == "__main__":
    main()