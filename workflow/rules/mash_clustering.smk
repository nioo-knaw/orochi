rule mash_sketch:
    input:
        forward = lambda wildcards: samples.loc[samples["sample"] == wildcards.sample].fq1.item(),
        rev = lambda wildcards: samples.loc[samples["sample"] == wildcards.sample].fq2.item(),
    output:
        temp(f"{outdir}/results/clustering/{{sample}}.msh")
    params:
       prefix = lambda wildcards: samples.loc[samples["sample"] == wildcards.sample, "sample"].item(),
       out_dir = f"{outdir}/results/clustering/"
    conda: "../envs/mash.yaml"
    shell: "mash sketch -k 27 -s 10000 -o {params.out_dir}{params.prefix} -r {input.forward} {input.rev}"

rule mash_paste:
    input:
        expand(f"{outdir}/results/clustering/{{sample}}.msh", sample=samples["sample"])
    output:
        temp("results/clustering/all.msh")
    conda: "../envs/mash.yaml"
    shell: "mash paste {rules.mash_sketch.params.out_dir}all {input}"

rule mash_dist:
    input: rules.mash_paste.output
    output: f"{outdir}/results/clustering/mashDistances.txt",
    conda: "../envs/mash.yaml"
    shell: "mash dist -t {input} {input} > {output}"

rule clustering:
    input: rules.mash_dist.output
    output: f"{outdir}/.samplesFileModified.checkpoint"
    params:
        clusterNumber = config["clustering"]
    conda: "../envs/python_clustering.yaml"
    shell:
        """
        python3 scripts/clustering.py -k {params.clusterNumber} && \
        touch {output}
        """