rule metabat:
    input: 
        contigs="{project}/assembly/{assembler}/{treatment}/{kmers}/assembly.fa",
        bam=expand("{{project}}/bamm/{{assembler}}/{{treatment}}/{{kmers}}/assembly.{sample}_forward_paired.bam", sample=config["data"]) 
    output:
        depth="{project}/binning/{assembler}/{treatment}/{kmers}/metabat/depth.txt",
        bin="{project}/binning/{assembler}/{treatment}/{kmers}/metabat/bin.1.fa"
    params:
        prefix="{project}/metabat/bin"
    log: "{project}/metabat/metabat.log"
    threads: 16
    run: 
        shell("/data/tools/metabat/dev/bin/jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}")
        shell("/data/tools/metabat/dev/bin/metabat -i {input.contigs} -a {output.depth} -o {params.prefix} --very-sensitive --numThreads {threads} --minContigByCorr 1500 --saveTNF saved.tnf --saveDistance saved.dist --verbose > {log}")

rule mmgenome_metabat:
    input:
        bin="{project}/metabat/bin.1.fa"
    output:
        "{project}/mmgenome/metabat.bins.txt"
    params:
        dir="{project}/metabat/"
    run:
        shell("echo -e 'scaffold\tbin' > {output}")
        shell("for file in `ls {params.dir}/*.fa` ; do noext=${{file%.fa}}; bin=$(basename ${{noext}}); awk -v bin=$bin '/^>/ {{split(substr($0,2),a,\":\"); print a[1] \"\\t\" bin;}}' $file;  done >> {output}")
        shell("sed --in-place -e s/[.]//g {output}")

rule checkm_lineage_metabat:
    input:
        bin="{project}/metabat/bin.1.fa"
    output:
        log="{project}/metabat/checkm/log.txt"
    params:
        indir="{project}/metabat/",
        outdir="{project}/metabat/checkm/"
    threads: 16
    shell: "set +u; source ~/.virtualenvs/groopm/bin/activate; source /data/tools/CheckM/0.9.7/env.sh; source /data/tools/pplacer/1.1/env.sh; set -u; python /data/tools/CheckM/0.9.7/bin/checkm lineage_wf -x fa -t {threads} -f {output.log} {params.indir} {params.outdir}"

