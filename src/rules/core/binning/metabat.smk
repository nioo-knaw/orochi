rule metabat:
    input: 
        contigs="scratch/assembly/megahit/minimus2/secondary.contigs.fasta",
        bam=expand("scratch/coverm/bamfiles/secondary.contigs.fasta.{sample}_R1.fastq.bam", sample=config["data"]) 
    output:
        depth="results/binning/metabat/metabat/depth.txt",
        bin="results/binning/metabat/metabat/bin.1.fa"
    params:
        prefix="results/binning/metabat/bin"
    log: "results/binning/metabat/metabat.log"
    threads: 16
    run: 
        shell("/data/tools/metabat/dev/bin/jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}")
        shell("/data/tools/metabat/dev/bin/metabat -i {input.contigs} -a {output.depth} -o {params.prefix} --very-sensitive --numThreads {threads} --minContigByCorr 1500 --saveTNF saved.tnf --saveDistance saved.dist --verbose > {log}")

rule mmgenome_metabat:
    input:
        bin="results/binning/metabat/bin.1.fa"
    output:
        "results/binning/mmgenome/metabat.bins.txt"
    params:
        dir="results/binning/metabat/"
    run:
        shell("echo -e 'scaffold\tbin' > {output}")
        shell("for file in `ls {params.dir}/*.fa` ; do noext=${{file%.fa}}; bin=$(basename ${{noext}}); awk -v bin=$bin '/^>/ {{split(substr($0,2),a,\":\"); print a[1] \"\\t\" bin;}}' $file;  done >> {output}")
        shell("sed --in-place -e s/[.]//g {output}")

rule checkm_lineage_metabat:
    input:
        bin="results/binning/metabat/bin.1.fa"
    output:
        log="results/binning/metabat/checkm/log.txt"
    params:
        indir="results/binning/metabat/",
        outdir="results/binning/metabat/checkm/"
    threads: 16
    shell: "set +u; source ~/.virtualenvs/groopm/bin/activate; source /data/tools/CheckM/0.9.7/env.sh; source /data/tools/pplacer/1.1/env.sh; set -u; python /data/tools/CheckM/0.9.7/bin/checkm lineage_wf -x fa -t {threads} -f {output.log} {params.indir} {params.outdir}"

