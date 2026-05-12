""" The rules related to Biosynthetic Gene Cluster (BGC) prediction and related analyses"""
from pathlib import Path
import json

def fungismash_input(wildcards):
    whokaryote_ckpt = checkpoints.whokaryote.get(sample=wildcards.sample_pool)
    euk_fasta = whokaryote_ckpt.output.euk_fasta

    if not fasta_has_records(euk_fasta):
        return {
            "contigs": euk_fasta,
        }

    augustify_ckpt = checkpoints.augustify.get(sample=wildcards.sample_pool)
    gff = augustify_ckpt.output.gff

    return {
        "db_ready": antismash_db_input(wildcards),
        "gff": gff,
        "contigs": euk_fasta,
    }

minsize_antismash = config['min_contig_antismash']
antismash_db = config["antismash_db"]
download_antismash_db = config.get("download_antismash_db", False)


def antismash_db_input(wildcards):
    if download_antismash_db:
        return os.path.join(antismash_db, ".download_complete")
    return []

rule download_antismash_databases:
    output:
        touch(os.path.join(antismash_db, ".download_complete"))
    conda:
        "../envs/antismash.yaml"
    params:
        dbdir=antismash_db,
    shell:
        """
        mkdir -p {params.dbdir}
        download-antismash-databases --database-dir {params.dbdir}
        """


rule antismash:
    input:
        db_ready=antismash_db_input,
        gff=f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_prokaryote_{minsize_antismash}.gff",
        contigs=f"{outdir}/results/04_gene_prediction/whokaryote/{{sample_pool}}/prokaryotes.fasta"
    output:
        html=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/index.html",
        json=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/bacterial.json"
    conda:
        "../envs/antismash.yaml"
    params:
        outdir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial",
        threads=config['threads'],
        database_dir=antismash_db

    shell:
        """
        test -d {params.database_dir} || \
            (echo "ERROR: antiSMASH database directory not found: {params.database_dir}" >&2; exit 1)
        
        antismash {input.contigs} \
        -c {params.threads} \
        --genefinding-gff3 {input.gff} \
        --output-dir {params.outdir} \
        --taxon bacteria \
        --output-basename bacterial \
        --cc-mibig \
        --cb-general \
        --cb-knownclusters \
        --databases {params.database_dir}
        """

rule fungismash:
    input:
        unpack(fungismash_input)

    output:
        html=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/index.html",
        json=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/fungal.json"

    log:
        f"{outdir}/logs/antismash/fungal/{{sample_pool}}.log"

    conda:
        "../envs/antismash.yaml"

    params:
        outdir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/",
        threads=config['threads'],
        database_dir=antismash_db

    run:
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


rule summarize_antismash_bacterial:
    input:
        json=rules.antismash.output.json
    output:
        tsv=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/bacterial_summary.tsv"
    conda:
        "../envs/antismash.yaml"
    params:
        antismash_dir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial",
        taxon="bacteria"
    script:
        "../scripts/summarize_antismash.py"


rule summarize_antismash_fungal:
    input:
        json=rules.fungismash.output.json
    output:
        tsv=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/fungal_summary.tsv"
    conda:
        "../envs/antismash.yaml"
    params:
        antismash_dir=f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal",
        taxon="fungi"
    script:
        "../scripts/summarize_antismash.py"


rule bigscape:
    input:
        "path/to/antismash_output"

    output:
        "path/to/output"

    shell:
        "bigscape -options"

rule itol_bgc:
    input:
        "path/to/input"
    output:
        "path/to/output"
