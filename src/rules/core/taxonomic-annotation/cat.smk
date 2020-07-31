rule CAT:
    input:
        "scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output:
        "scratch/annotation/CAT/secondary.contigs.ORF2LCA.txt",
        "scratch/annotation/CAT/secondary.contigs.contig2classification.txt",
    params:
        db=config['CAT_database'],
        tax=config['CAT_taxonomy'],
        tmp=config['tmpdir'],
        prefix="scratch/annotation/CAT/secondary.contigs"

    log: "scratch/annotation/CAT/secondary.contigs.log"
    conda:
        "../../../envs/cat.yaml"
    threads: 16
    shell:
        """CAT contigs -c {input} -o {params.prefix} -d {params.db} -t {params.tax} --nproc {threads} --sensitive --tmpdir {params.tmp} --no_log > {log} 2>&1"""


rule CAT_add_names:
    input:
        "scratch/annotation/CAT/secondary.contigs.ORF2LCA.txt",
    output:
        "scratch/annotation/CAT/secondary.contigs.classification.txt",
    params:
        tax=config['CAT_taxonomy'],
    log: "scratch/annotation/CAT/secondary.contigs.names.log"
    conda:
        "../../../envs/cat.yaml"
    shell: "CAT add_names -i {input} -o {output} -t {params.tax} --only_official > {log} 2>&1"

rule CAT_summarize:
    input:
        contigs = "scratch/assembly/megahit/minimus2/secondary.contigs.fasta",
        classification = "scratch/annotation/CAT/secondary.contigs.classification.txt"
    output:
        "scratch/annotation/CAT/secondary.contigs.classification.summary.txt"
    log: "scratch/annotation/CAT/secondary.contigs.summary.log"
    conda:
        "../../../envs/cat.yaml"
    shell: "CAT summarise -c {input.contigs} -i {input.classification} -o {output} > {log} 2>&1"
