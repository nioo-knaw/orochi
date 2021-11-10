rule pre_metabat:
    input: "scratch/assembly/megahit/minimus2/all.merged.contigs.fasta"
    output: "scratch/assembly/megahit/minimus2/filtered_assembly.fasta"
    log: "logs/binning/metabat/pre_metabat.log"
    conda: "../../../envs/orochi-base.yaml"
    threads: 16
    shell: "python src/scripts/assembly_rmdup.py"

rule metabat:
    input: 
        contigs="scratch/assembly/megahit/minimus2/filtered_assembly.fasta",
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
        mkdir results/binning/metabat
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        metabat2 -i {input.contigs} -a {output.depth} -o {params.prefix} -v > {log}
        """

rule mmgenome_metabat:
    input:
        bin="results/binning/metabat/bins/bin.1.fa"
    output:
        "results/binning/metabat/metabat.bins.txt"
    params:
        dir="results/binning/metabat/bins"
    run:
        shell("echo -e 'scaffold\tbin' > {output}")
        shell("for file in `ls {params.dir}/*.fa` ; do noext=${file%.fa}; bin=$(basename ${noext}); awk -v bin=$bin '/^>/ {split(substr($0,2),a; print a[1] "\t" bin;}' $file;  done >> {output}")
        shell("sed --in-place -e s/[.]//g {output}")
