rule coverage:
    input:
        forward = expand("scratch/host_filtering/{sample}_R1.fastq", sample=config["data"]),
        reverse = expand("scratch/host_filtering/{sample}_R2.fastq", sample=config["data"]),
        assembly=expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        "results/stats/coverage.tsv"
    conda:
        "../../../envs/coverm.yaml"
    threads: 16
    log:
        "results/logs/coverm/coverage.log"
    shell:
        "coverm contig --mapper bwa-mem --methods mean --reference {input.assembly} -1 {input.forward} -2 {input.reverse} -t {threads} 2> {log} > {output}"

#TO DO: Add stderr log?
