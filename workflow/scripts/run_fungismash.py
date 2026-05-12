from pathlib import Path
import json

outdir_path = Path(params.outdir)
outdir_path.mkdir(parents=True, exist_ok=True)

log_path = Path(log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w") as log_handle:
    log_handle.write(f"[fungismash] sample_pool={wildcards.sample_pool}\n")
    log_handle.write(f"[fungismash] input={list(input)}\n")
    log_handle.write(f"[fungismash] output.html={output.html}\n")
    log_handle.write(f"[fungismash] output.json={output.json}\n")
    log_handle.write(f"[fungismash] outdir={params.outdir}\n")
    log_handle.write(f"[fungismash] database_dir={params.database_dir}\n")

# Case 1: no eukaryotic sequences
if not fasta_has_records(input.contigs):
    reason = "No eukaryotic sequences are detected"

    with open(log_path, "a") as log_handle:
        log_handle.write(f"[fungismash] {reason}. Writing placeholder outputs.\n")

    with open(output.html, "w") as handle:
        handle.write(
            "<html>\n"
            "<head><title>No eukaryotic sequences</title></head>\n"
            "<body>\n"
            "<h1>No eukaryotic sequences are detected</h1>\n"
            "</body>\n"
            "</html>\n"
        )

    with open(output.json, "w") as handle:
        json.dump(
            {
                "sample": wildcards.sample_pool,
                "taxon": "fungi",
                "status": "skipped",
                "reason": reason,
                "records": [],
            },
            handle,
            indent=2,
        )

# Case 2: eukaryotic contigs exist, but augustify produced no gene annotation
elif "gff" in input.keys() and not gff_has_features(input.gff):
    reason = "Eukaryotic sequences are detected, but no fungal gene annotations are available"

    with open(log_path, "a") as log_handle:
        log_handle.write(f"[fungismash] {reason}. Writing placeholder outputs.\n")
        log_handle.write(f"[fungismash] gff={input.gff}\n")

    with open(output.html, "w") as handle:
        handle.write(
            "<html>\n"
            "<head><title>No fungal gene annotations</title></head>\n"
            "<body>\n"
            "<h1>Eukaryotic sequences are detected, but no fungal gene annotations are available</h1>\n"
            "<p>fungiSMASH was skipped because the augustify GFF file contains no gene, mRNA, or CDS features.</p>\n"
            "</body>\n"
            "</html>\n"
        )

    with open(output.json, "w") as handle:
        json.dump(
            {
                "sample": wildcards.sample_pool,
                "taxon": "fungi",
                "status": "skipped",
                "reason": reason,
                "records": [],
            },
            handle,
            indent=2,
        )

# Case 3: normal fungiSMASH run
else:
    shell(
        r"""
        echo "[fungismash] Eukaryotic sequences and gene annotations detected. Running antiSMASH." >> {log} 2>&1

        test -d {params.database_dir} || \
            (echo "ERROR: antiSMASH database directory not found: {params.database_dir}" >> {log} 2>&1; exit 1)

        echo "[fungismash] FASTA record count:" >> {log} 2>&1
        grep -c "^>" {input.contigs} >> {log} 2>&1 || true

        echo "[fungismash] GFF feature preview:" >> {log} 2>&1
        awk '$0 !~ /^#/ {{print; count++}} count==5 {{exit}}' {input.gff} >> {log} 2>&1 || true

        antismash {input.contigs} \
            -c {params.threads} \
            --genefinding-gff3 {input.gff} \
            --output-dir {params.outdir} \
            --taxon fungi \
            --cassis \
            --output-basename fungal \
            --cc-mibig \
            --cb-general \
            --cb-knownclusters \
            --databases {params.database_dir} \
            >> {log} 2>&1

        echo "[fungismash] antiSMASH finished." >> {log} 2>&1

        test -s {output.html} || \
            (echo "ERROR: expected output missing or empty: {output.html}" >> {log} 2>&1; exit 1)

        test -s {output.json} || \
            (echo "ERROR: expected output missing or empty: {output.json}" >> {log} 2>&1; exit 1)
        """
    )