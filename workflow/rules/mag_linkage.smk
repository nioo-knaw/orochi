""" Rules related to reconstructing 16S rRNA gene sequences and linking them to MAGs"""

rule phyloflash:
    input:
        forward=f"{'outdir'}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
        reverse=f"{'outdir'}/results/03_assembly/coassembly/pools/{{sample_pool}}_rev.fastq.gz"

    output:
        phyloflash_dir=f"{'outdir'}/results/07_phyloflash/{{sample_pool}}",
        phyloflash_report=f"{'outdir'}/results/07_phyloflash/{{sample_pool}}/{{sample_pool}}_phyloFlash.report.csv",
        phyloflash_fasta=f"{'outdir'}/results/07_phyloflash/{{sample_pool}}/{{sample_pool}}.all.final.fasta"

    conda:
        "../envs/phyloflash.yaml"

    params:
        db=config["phyloflash_db"],
        threads=config["threads"]

    shell:
        "phyloFlash.pl -o {output.phyloflash_dir} -dbhome {params.db} -lib {wildcards.sample_pool} \
         -CPUs {params.threads} -read1 {input.forward} -read2 {input.reverse}"