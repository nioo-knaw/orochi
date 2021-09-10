rule metabat:
    input: 
        contigs="scratch/assembly/megahit/minimus2/all.merged.contigs.fasta",
        bam=expand("scratch/coverm/bamfiles/all.merged.contigs.fasta.{sample}_R1.fastq.bam", sample=config["data"]) 
    output:
        depth="results/binning/metabat/depth.txt",
        bin="results/binning/metabat/bins/bin.1.fa"
    params:
        prefix="results/binning/metabat/bins/bin"
    log: "logs/binning/metabat/metabat.log"
    conda:
        "../../../envs/metabat.yaml"
    threads: 16
    shell: 
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        metabat -i {input.contigs} -a {output.depth} -o {params.prefix} --minContig 1500 -v > {log}
        """

rule mmgenome_metabat:
    input:
        bin="results/binning/metabat/bins/bin.1.fa"
    output:
        "results/binning/metabat/metabat.bins.txt"
    params:
        dir="results/binning/metabat/bins"
    log: "logs/metabat/mmgenome_metabat.log"
    run:
        shell("echo -e 'scaffold\tbin' > {output}")
        shell("for file in `ls {params.dir}/*.fa` ; do noext=${{file%.fa}}; bin=$(basename ${{noext}}); awk -v bin=$bin '/^>/ {{split(substr($0,2),a,\":\"); print a[1] \"\\t\" bin;}}' $file;  done >> {output}")
        shell("sed --in-place -e s/[.]//g {output} 2> {log}")
