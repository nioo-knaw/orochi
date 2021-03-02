rule checkm_lineage:
    input: "results/binning/{binner}/"
    output:
        "results/binning/{binner}/checkm/completeness.tsv",
        "results/binning/{binner}/checkm/concatenated.fasta"
    params:
        outdir="results/binning/{binner}/checkm"
    threads: 40
    conda: "../../../envs/checkm.yaml"
    shell: "checkm lineage_wf -t {threads} -x fa {input} {params.outdir}"
