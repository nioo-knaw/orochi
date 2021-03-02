rule groopm:
    input:
        contigs="scratch/assembly/megahit/minimus2/secondary.contigs.fasta",
        bam=expand("scratch/coverm/bamfiles/secondary.contigs.fasta.{sample}_R1.fastq.bam", sample=config["data"])
    output: "results/binning/groopm/bin1.fasta"
    params:
        outdir="results/binning/groopm/"
    conda: "../../../envs/groomp.yaml"
    shell:
        """
        groopm parse db.gm {input.contigs} scratch/coverm/bamfiles/*.bam
        groopm core db.gm
        groopm refine db.gm
        groopm recruit db.gm
        groopm extract db.gm {input.contigs} -o {params.outdir}
        """
