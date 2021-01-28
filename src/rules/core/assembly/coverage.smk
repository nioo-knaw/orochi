rule coverage:
    input:
        forward = expand("scratch/treatment/{treatment}_forward.fastq", treatment=config["treatment"]),
        reverse = expand("scratch/treatment/{treatment}_reverse.fastq", treatment=config["treatment"]),
        #assembly=expand("scratch/assembly/megahit/{treatment}/{kmers}/assembly.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
        assembly="scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output:
        "results/stats/coverage.tsv"
        #"scratch/coverm/coverage.tsv"
    params:
        bamdir="scratch/coverm/bamfiles"
    conda:
        "../../../envs/coverm.yaml"
    threads: 16
    shell:
        "coverm contig --mapper bwa-mem --methods mean -c {input.forward} {input.reverse} --reference {input.assembly} --bam-file-cache-directory {params.bamdir} -t {threads} -o {output}"

#TO DO: Add stderr log?
