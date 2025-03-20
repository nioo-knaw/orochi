""" Rules related to reconstructing 16S rRNA gene sequences and linking them to MAGs"""

rule phyloflash:
    input:
        forward_reads=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
        reverse_reads=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_rev.fastq.gz"

    output:
        phyloflash_out=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.phyloFlash.tar.gz",
        # phyloflash_report=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}_phyloFlash.report.csv",
        # phyloflash_fasta=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.all.final.fasta"

    conda:
        "../envs/phyloflash.yaml"

    params:
        db=config["phyloflash_db"],
        threads=config["threads"],
        phylo_dir=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/"

    shell: # We have to zip the phyloflash output and move it becaues it will be stored in the WD otherwise.
        "phyloFlash.pl -dbhome {params.db} -lib {wildcards.sample_pool} -zip \
         -CPUs {params.threads} -read1 {input.forward_reads} -read2 {input.reverse_reads}; \ mv {wildcards.sample_pool}.phyloFlash.* {params.phylo_dir}"

rule unzip_phyloflash:
    input:
        phyloflash_tar=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.phyloFlash.tar.gz"
    output:
        phyloflash_output=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.all.final.fasta"
    shell:
        "tar -xzf {input.phyloflash_tar}"

rule unzip_reads:
    input:
        forward_reads=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
        reverse_reads=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_rev.fastq.gz"
    output:
        forward_out=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_forward.fastq",
        reverse_out=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_reverse.fastq"
    shell:
        "gunzip -c {input.forward_reads} > {output.forward_out}; gunzip -c {input.reverse_reads} > {output.reverse_out}"

rule rename_reads:
    input:
        forward_reads=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_forward.fastq",
        reverse_reads=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_reverse.fastq"
    output:
        forward_renamed=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R1.fastq",
        reverse_renamed=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R2.fastq"
    conda:
        "../envs/markerMAG.yaml"
    params:
        threads=config["threads"]
    shell:
        "MarkerMAG rename_reads -r1 {input.forward_reads} -r2 {input.reverse_reads} -p {wildcards.sample_pool} -fq \
        -t {params.threads}"

rule fastq_2_fasta:
    input:
        forward_reads=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R1.fastq",
        reverse_reads=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R2.fastq"
    output:
        fasta_forward=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R1.fasta",
        fasta_reverse=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R2.fasta"
    conda:
        "../envs/seqkit.yaml"
    params:
        threads=config["threads"]
    shell:
        "seqkit fq2fa {input.forward_reads} -o {input.reverse_reads} --threads {params.threads}"


rule markermag_link:
    input:
        forward_reads=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R1.fasta",
        reverse_reads=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R2.fasta",
        phyloflash=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.all.final.fasta",
        mag_fasta=directory(f"{outdir}/results/06_binning/drep/dereplicated_genomes")
    output:
        markerMAG_link=f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag/{{sample_pool}}_markerMAG_link.csv"
    conda:
        "../envs/markerMAG.yaml"
    params:
        threads=config["threads"]
    shell:
        "MarkerMAG link -p {wildcards.sample_pool} -r1 {input.forward_reads} -r2 {input.reverse_reads} \
        -marker {input.phyloflash} -mag {input.mag_fasta} -o {output.markerMAG_link} -x fa -t {params.threads}"