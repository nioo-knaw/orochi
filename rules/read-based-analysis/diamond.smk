rule diamond_fraggenescan:
    input:
        fasta="{project}/read-based/fraggenescan/{sample}_forward_paired.ffn"
    output:
        tsv="{project}/read-based/diamond/{sample}_forward_paired.diamond.nr.daa",
    params:
        reference=config["diamond_database"],
        version="0.9.14",
        output="{project}/read-based/diamond/{sample}_forward_paired.diamond.nr",
        format="tab",
        tmp="/tmp",
        megan_version=config['megan_version'],
        megan_mapping=config['megan_mapping']
    priority: 20
    threads: 32
    shell: "/data/tools/diamond/{params.version}/bin/diamond blastx --sensitive -c 1 -d {params.reference} -t {params.tmp} -p {threads} -q {input} -a {params.output}"
