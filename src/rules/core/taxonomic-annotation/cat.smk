rule CAT:
    input:
        # TODO: Decide what is the input here
        expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        temporary("scratch/annotation/CAT/assembly.alignment.diamond.gz"),
        "scratch/annotation/CAT/assembly.predicted_proteins.gff",
        "scratch/annotation/CAT/assembly.predicted_proteins.faa",
        "scratch/annotation/CAT/assembly.ORF2LCA.txt",
        "scratch/annotation/CAT/assembly.contig2classification.txt",

    params:
        db=config['CAT_database'],
        tax=config['CAT_taxonomy'],
        tmp=config['tmpdir'],
        prefix="scratch/annotation/CAT/assembly"

    log: "scratch/annotation/CAT/assembly.log"
    conda:
        "../../../envs/cat.yaml"
    threads: 16
    shell:
        """CAT contigs -c {input} -o {params.prefix} -d {params.db} -t {params.tax} --nproc {threads} --sensitive --tmpdir {params.tmp} --no_log > {log} 2>&1"""


rule CAT_add_names:
    input:
        "scratch/annotation/CAT/assembly.contig2classification.txt",
    output:
        "scratch/annotation/CAT/assembly.classification.txt",
    params:
        tax=config['CAT_taxonomy'],
    log: "scratch/annotation/CAT/assembly.names.log"
    conda:
        "../../../envs/cat.yaml"
    shell: "CAT add_names -i {input} -o {output} -t {params.tax} --only_official > {log} 2>&1"

rule CAT_summarize:
    input:
        contigs = expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"]),
        classification = "scratch/annotation/CAT/assembly.classification.txt"
    output:
        "results/annotation/CAT/assembly.classification.summary.txt"
    log: "scratch/annotation/CAT/assembly.summary.log"
    conda:
        "../../../envs/cat.yaml"
    shell: "CAT summarise -c {input.contigs} -i {input.classification} -o {output} > {log} 2>&1"
