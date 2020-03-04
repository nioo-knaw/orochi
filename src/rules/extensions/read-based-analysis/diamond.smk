rule diamond_fraggenescan:
    input:
        fasta="scratch/read-based/fraggenescan/{sample}_forward_paired.ffn"
    output:
        tsv="scratch/read-based/diamond/{sample}_forward_paired.diamond.nr.daa",
    params:
        reference=config["diamond_database"],
        version="0.9.14",
        output="scratch/read-based/diamond/{sample}_forward_paired.diamond.nr",
        format="tab",
        tmp="/tmp",
        megan_version=config['megan_version'],
        megan_mapping=config['megan_mapping']
    priority: 20
    threads: 32
    shell: "/data/tools/diamond/{params.version}/bin/diamond blastx -c 1 -d {params.reference} -t {params.tmp} -p {threads} -q {input} -a {params.output}"

rule diamond_taxonomy_and_kegg:
    input:
        "scratch/read-based/diamond/{sample}_forward_paired.diamond.nr.daa"
    output:
        taxonomy="scratch/read-based/diamond/{sample}_forward_paired.diamond.nr-taxonomy.txt",
        kegg="scratch/read-based/diamond/{sample}_forward_paired.diamond.nr-kegg.txt"
    threads: 32
    shell: "/data/tools/megan-ue/6.10.8/tools/blast2lca -i {input} -f DAA -ms 50 -me 0.01 -top 50 -a2t /data/db/megan/prot_acc2tax-Mar2018X1.abin -a2kegg /data/db/megan/acc2kegg-Dec2017X1-ue.abin --kegg"

