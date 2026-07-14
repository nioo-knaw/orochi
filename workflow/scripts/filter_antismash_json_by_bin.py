#!/usr/bin/env python3
"""Filter an antiSMASH results JSON down to the records belonging to one bin.

antiSMASH's results JSON (e.g. bacterial.json) stores one entry per contig
in its "records" list, each with an "id" (and sometimes a separate
"original_id", when antiSMASH had to rewrite the id, e.g. for GenBank's
16-character LOCUS limit). Keeping only the records whose id/original_id
belongs to one bin gives a JSON that antiSMASH's own `--reuse-results` mode
can turn back into a full HTML report for just that bin, without re-running
any analysis modules.
"""

import json
import sys
from pathlib import Path


def load_contig_to_bin(contig2bin_file):
    """Parse a DAS_Tool-style contig2bin TSV (no header: contig\tbin)."""
    contig_to_bin = {}
    with open(contig2bin_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            contig_to_bin[parts[0]] = parts[1]
    return contig_to_bin


def filter_json_by_bin(antismash_json, contig2bin_file, bin_id, output_json):
    """Write a copy of antismash_json containing only records for bin_id.

    Returns the number of records kept.
    """
    contig_to_bin = load_contig_to_bin(contig2bin_file)
    bin_contigs = {c for c, b in contig_to_bin.items() if b == bin_id}

    if not bin_contigs:
        raise ValueError(f"No contigs assigned to bin '{bin_id}' in {contig2bin_file}")

    with open(antismash_json) as fh:
        data = json.load(fh)

    kept = []
    for record in data["records"]:
        record_id = record.get("id")
        original_id = record.get("original_id", record_id)
        if record_id in bin_contigs or original_id in bin_contigs:
            kept.append(record)

    if not kept:
        raise ValueError(
            f"Bin '{bin_id}' has {len(bin_contigs)} contig(s) assigned in "
            f"{contig2bin_file}, but none of them matched a record id/"
            f"original_id in {antismash_json}. Check that contig naming is "
            "consistent between the assembly used for binning and the one "
            "fed to antiSMASH."
        )

    filtered = dict(data)
    filtered["records"] = kept
    filtered["input_file"] = f"{bin_id}_reuse"

    if isinstance(data.get("timings"), dict):
        kept_ids = {r["id"] for r in kept}
        filtered["timings"] = {k: v for k, v in data["timings"].items() if k in kept_ids}

    output_json = Path(output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with output_json.open("w") as fh:
        json.dump(filtered, fh)

    return len(kept)


if __name__ == "__main__":
    log = open(snakemake.log[0], "w")
    sys.stdout = log
    sys.stderr = log

    n_records = filter_json_by_bin(
        antismash_json=snakemake.input.antismash_json,
        contig2bin_file=snakemake.input.contig2bin,
        bin_id=snakemake.params.bin_id,
        output_json=snakemake.output.filtered_json,
    )

    print(f"Kept {n_records} record(s) for bin {snakemake.params.bin_id}")