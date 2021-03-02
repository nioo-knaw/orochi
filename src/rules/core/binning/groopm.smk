rule groopm:
    input:
        contigs="scratch/assembly/megahit/minimus2/secondary.contigs.fasta",
        bam=expand("scratch/coverm/bamfiles/secondary.contigs.fasta.{sample}_R1.fastq.bam", sample=config["data"])
    output: "results/binning/groopm/bin1.fasta"
    params:
        db="results/binning/groopm/db.gm",
        outdir="results/binning/groopm/"
    conda: "../../../envs/groopm.yaml"
    threads: 40
    shell:
        """
        groopm parse {params.db} {input.contigs} scratch/coverm/bamfiles/*.bam -t {threads}
        groopm core {params.db} -t {threads}
        groopm refine {params.db} -t {threads}
        groopm recruit {params.db} -t {threads}
        groopm extract {params.db} {input.contigs} -o {params.outdir} -t {threads}
        """
