rule checkm_lineage_metabat:
    input:
        "results/binning/metabat/bins/bin.1.fa"
    output:
        "results/binning/metabat/checkm/completeness.tsv",
        "results/binning/metabat/checkm/concatenated.fasta"
    params:
        indir=lambda wildcards, input: os.path.dirname(os.path.dirname(str(input))),
        outdir=lambda wildcards, output: os.path.dirname(output[0])
    log:
        "logs/binning/checkm_lineage_metabat.log"
    threads: 40
    conda: "../../../envs/checkm.yaml"
    shell: "checkm lineage_wf -t {threads} -x fa {params.indir} {params.outdir} &> {log}"
