"""
rule checkm_lineage:
    input:
        "results/binning/{binner}/bin1.fasta",
        "results/binning/{binner}/bin.1.fa"
    output:
        "results/binning/{binner}/checkm/completeness.tsv",
        "results/binning/{binner}/checkm/concatenated.fasta"
    params:
        indir="results/binning/{binner}",
        outdir="results/binning/{binner}/checkm"
    threads: 40
    conda: "../../../envs/checkm.yaml"
    shell: "checkm lineage_wf -t {threads} -x fa {params.indir} {params.outdir}"
"""
