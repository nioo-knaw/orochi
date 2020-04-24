rule antismash:
    input:
        "scratch/assembly/megahit/secondary.contigs.fasta"
    output:
        "scratch/annotation/antismash/secondary.contigs.gbk"
    params:
        outdir="scratch/annotation/antismash/"
    conda:
        "../../../envs/antismash.yaml"
    shell:
        "antismash --genefinding-tool prodigal-m --output-dir {params.outdir} --cpus {threads} {input}"
