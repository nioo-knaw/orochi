rule antismash:
    input:
        "scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output:
        "scratch/annotation/antismash/secondary.contigs.gbk"
    params:
        outdir="scratch/annotation/antismash/"
    conda:
        "../../../envs/antismash.yaml"
    shell:
        "antismash --cb-general --cb-knownclusters --cb-subclusters --asf --genefinding-tool prodigal-m --output-dir {params.outdir} --cpus {threads} {input}"
