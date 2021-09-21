rule metabat:
    input: 
        contigs="scratch/assembly/megahit/minimus2/all.merged.contigs.fasta",
        bam=expand("scratch/coverm/bamfiles/all.merged.contigs.fasta.{sample}_R1.fastq.bam", sample=config["data"]) 
    output:
        depth="results/binning/metabat/depth.txt",
        bin="results/binning/metabat/bins/bin.1.fa"
    params:
        prefix=lambda wildcards, output: os.path.join(os.path.dirname(output[1]), "bin")
    log: "logs/binning/metabat/metabat.log"
    conda:
        "../../../envs/metabat.yaml"
    threads: 16
    shell: 
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        metabat -i {input.contigs} -a {output.depth} -o {params.prefix} --minContig 1500 -v > {log}
        """

