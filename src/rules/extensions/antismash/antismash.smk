rule antismash:
    input:
        "scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output:
        "scratch/annotation/antismash/secondary.contigs.gbk",
        "scratch/annotation/antismash/secondary.contigs.json"
    params:
        outdir="scratch/annotation/antismash/"
    conda:
        "../../../envs/antismash.yaml"
    threads: 16
    shell:
        "antismash --cb-general --cb-knownclusters --cb-subclusters --asf --genefinding-tool prodigal-m --output-dir {params.outdir} --cpus {threads} {input}"


rule get_bgcs:
    input:
        "scratch/annotation/antismash/secondary.contigs.json"
    output:
        "scratch/annotation/antismash/bgcs.fasta"
    script:
        "../../../scripts/antismash_get_bgcs.py"
