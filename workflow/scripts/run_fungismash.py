#!/usr/bin/env python3

from pathlib import Path
import json
import subprocess
import sys
import shutil

def fasta_has_records(fasta):
    """
    Return True if a FASTA file exists and contains at least one sequence header.
    """
    fasta = Path(fasta)

    if not fasta.exists() or fasta.stat().st_size == 0:
        return False

    with open(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                return True

    return False


def gff_has_features(gff):
    """
    Return True if the GFF contains at least one gene, mRNA, or CDS feature.
    """
    gff = Path(gff)

    if not gff.exists() or gff.stat().st_size == 0:
        return False

    with open(gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2] in {"gene", "mRNA", "CDS"}:
                return True

    return False


def write_placeholder(html, json_out, sample, title, reason, detail=None):
    """
    Write placeholder fungiSMASH HTML and JSON outputs.
    """
    html = Path(html)
    json_out = Path(json_out)

    html.parent.mkdir(parents=True, exist_ok=True)
    json_out.parent.mkdir(parents=True, exist_ok=True)

    detail_html = ""
    if detail:
        detail_html = f"<p>{detail}</p>\n"

    with open(html, "w") as handle:
        handle.write(
            "<html>\n"
            f"<head><title>{title}</title></head>\n"
            "<body>\n"
            f"<h1>{reason}</h1>\n"
            f"{detail_html}"
            "</body>\n"
            "</html>\n"
        )

    with open(json_out, "w") as handle:
        json.dump(
            {
                "sample": sample,
                "taxon": "fungi",
                "status": "skipped",
                "reason": reason,
                "records": [],
            },
            handle,
            indent=2,
        )


def reset_output_dir(outdir):
    """
    Remove and recreate the antiSMASH output directory.
    """
    outdir_path = Path(outdir).resolve()

    if outdir_path.name != "fungal":
        raise ValueError(f"Refusing to remove unexpected output directory: {outdir_path}")

    if outdir_path.exists():
        shutil.rmtree(outdir_path)

    outdir_path.mkdir(parents=True, exist_ok=True)


def main():
    sample = snakemake.wildcards.sample_pool

    contigs = snakemake.input["contigs"]
    html = snakemake.output["html"]
    json_out = snakemake.output["json"]

    outdir = snakemake.params["outdir"]
    threads = str(snakemake.params["threads"])
    database_dir = snakemake.params["database_dir"]

    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)

    print(f"[fungismash] sample_pool={sample}")
    print(f"[fungismash] input={list(snakemake.input)}")
    print(f"[fungismash] output.html={html}")
    print(f"[fungismash] output.json={json_out}")
    print(f"[fungismash] outdir={outdir}")
    print(f"[fungismash] database_dir={database_dir}")

    # Case 1: no eukaryotic sequences.
    # Important: do this before accessing snakemake.input["gff"].
    if not fasta_has_records(contigs):
        reason = "No eukaryotic sequences are detected"
        print(f"[fungismash] {reason}. Writing placeholder outputs.")
        reset_output_dir(outdir)
        write_placeholder(
            html=html,
            json_out=json_out,
            sample=sample,
            title="No eukaryotic sequences",
            reason=reason,
        )

        return

    # Case 2a: eukaryotic contigs exist, but no GFF was passed.
    if "gff" not in snakemake.input.keys():
        reason = "Eukaryotic sequences are detected, but no fungal gene annotations are available"
        detail = "fungiSMASH was skipped because no augustify GFF file was provided."

        print(f"[fungismash] {reason}. Writing placeholder outputs.")
        reset_output_dir(outdir)
        write_placeholder(
            html=html,
            json_out=json_out,
            sample=sample,
            title="No fungal gene annotations",
            reason=reason,
            detail=detail,
        )

        return

    gff = snakemake.input["gff"]
    print(f"[fungismash] gff={gff}")

    # Case 2b: eukaryotic contigs exist, but augustify produced no gene annotation.
    if not gff_has_features(gff):
        reason = "Eukaryotic sequences are detected, but no fungal gene annotations are available"
        detail = (
            "fungiSMASH was skipped because the augustify GFF file contains "
            "no gene, mRNA, or CDS features."
        )

        print(f"[fungismash] {reason}. Writing placeholder outputs.")
        reset_output_dir(outdir)
        write_placeholder(
            html=html,
            json_out=json_out,
            sample=sample,
            title="No fungal gene annotations",
            reason=reason,
            detail=detail,
        )

        return

    # Case 3: normal fungiSMASH run.
    if not Path(database_dir).is_dir():
        raise FileNotFoundError(f"antiSMASH database directory not found: {database_dir}")

    cmd = [
        "antismash",
        str(contigs),
        "-c", str(threads),
        "--genefinding-gff3", str(gff),
        "--output-dir", str(outdir),
        "--taxon", "fungi",
        "--cassis",
        "--output-basename", "fungal",
        "--tfbs",
        "--cc-mibig",
        "--cb-general",
        "--cb-knownclusters",
        "--databases", str(database_dir),
    ]

    print("[fungismash] Eukaryotic sequences and gene annotations detected. Running antiSMASH.")
    print("[fungismash] command: " + " ".join(cmd))
    reset_output_dir(outdir)

    subprocess.run(
        cmd,
        stdout=sys.stdout,
        stderr=sys.stderr,
        check=True,
    )

    if not Path(html).exists() or Path(html).stat().st_size == 0:
        raise FileNotFoundError(f"Expected output missing or empty: {html}")

    if not Path(json_out).exists() or Path(json_out).stat().st_size == 0:
        raise FileNotFoundError(f"Expected output missing or empty: {json_out}")

    print("[fungismash] antiSMASH finished successfully.")


if __name__ == "__main__":
    main()