""" Rules related to reconstructing 16S rRNA gene sequences and linking them to MAGs"""

rule phyloflash:
    input:
        forward=f"{'outdir'}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
        reverse=f"{'outdir'}/results/03_assembly/coassembly/pools/{{sample_pool}}_rev.fastq.gz"

    output:
        phyloflash_dir=f"{'outdir'}/results/07_maglinkage/{{sample_pool}}/phyloflash",
        phyloflash_report=f"{'outdir'}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}_phyloFlash.report.csv",
        phyloflash_fasta=f"{'outdir'}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.all.final.fasta"

    conda:
        "../envs/phyloflash.yaml"

    params:
        db=config["phyloflash_db"],
        threads=config["threads"],
        outdir=f"{'outdir'}/results/07_maglinkage/{{sample_pool}}/phyloflash/"

    shell:
        "phyloFlash.pl -o {output.phyloflash_dir} -dbhome {params.db} -lib {wildcards.sample_pool} \
         -CPUs {params.threads} -read1 {input.forward} -read2 {input.reverse}; mv {wildcards.sample_pool}.* {params.outdir}"

rule blast_db:
    input:
        fasta=f"{'outdir'}/results/07_maglinkage/{{sample_pool}}/blastdb/{{sample_pool}}.all.final.fasta"

    output:
        db=f"{'outdir'}/results/07_maglinkage/{{sample_pool}}/blastdb/{{sample_pool}}.all.final.fasta.nhr"

    conda:
        "../envs/blast.yaml"

    params:
        outfile=f"{'outdir'}/results/07_maglinkage/{{wildcards.sample_pool}}/blastdb/{{wildcards.sample_pool}}.all.final"

    shell:
        "makeblastdb -in {input.fasta} -dbtype nucl -out {params.outfile}"