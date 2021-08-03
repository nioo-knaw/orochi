rule CAT:
    input:
        expand("scratch/assembly/{assembler}/minimus2/secondary.contigs.fasta", assembler=config["assembler"])
        #expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        temporary("results/annotation/CAT/assembly.alignment.diamond.gz"),
        "results/annotation/CAT/assembly.predicted_proteins.gff",
        "results/annotation/CAT/assembly.predicted_proteins.faa",
        "results/annotation/CAT/assembly.ORF2LCA.txt",
        "results/annotation/CAT/assembly.contig2classification.txt",
    params:
        db=config['CAT_database'],
        tax=config['CAT_taxonomy'],
        tmp=config['tmpdir'],
        prefix=lambda wildcards, output: output[0][:-3]
    log: "logs/annotation/CAT/assembly.log"
    conda:
        "../../../envs/cat.yaml"
    threads: 80
    shell:
        """CAT contigs -c {input} -o {params.prefix} -d {params.db} -t {params.tax} --nproc {threads} --sensitive --force --compress --verbose --index_chunks 1 --top 11 --I_know_what_Im_doing"""


rule CAT_add_names:
    input:
        "results/annotation/CAT/assembly.contig2classification.txt",
    output:
        "results/annotation/CAT/assembly.classification.txt",
    params:
        tax=config['CAT_taxonomy'],
    log: "logs/annotation/CAT/assembly.names.log"
    conda:
        "../../../envs/cat.yaml"
    shell: "CAT add_names -i {input} -o {output} -t {params.tax} --only_official > {log} 2>&1"

rule CAT_summarize:
    input:
        #contigs = expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"]),
        contigs = expand("scratch/assembly/{assembler}/minimus2/secondary.contigs.fasta", assembler=config["assembler"]),
        classification = "results/annotation/CAT/assembly.classification.txt"
    output:
        protected("results/annotation/CAT/assembly.classification.summary.txt")
    log: "logs/annotation/CAT/assembly.summary.log"
    conda:
        "../../../envs/cat.yaml"
    shell: "CAT summarise -c {input.contigs} -i {input.classification} -o {output} > {log} 2>&1"
