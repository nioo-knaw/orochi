rule coverage:
    input:
        forward = expand("scratch/host_filtering/{sample}_R1.fastq", sample=config["data"]),
        reverse = expand("scratch/host_filtering/{sample}_R2.fastq", sample=config["data"]),
        assembly=expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        "results/stats/coverage.tsv"
    params:
        bamdir="scratch/coverm/bamfiles"
    conda:
        "../../../envs/coverm.yaml"
    threads: 16
    shell:
        "coverm contig --mapper bwa-mem --methods mean -c {input.forward} {input.reverse} --reference {input.assembly} --bam-file-cache-directory {params.bamdir} -t {threads} -o {output}"

#TO DO: Add stderr log?
