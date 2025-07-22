""" Rules related to reconstructing 16S rRNA gene sequences and linking them to MAGs"""

rule phyloflash:
    input:
        forward_reads=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
        reverse_reads=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_rev.fastq.gz"

    output:
        phyloflash_out=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.phyloFlash.tar.gz",
        phyloflash_done=touch(f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/phyloflash.done")
        # phyloflash_report=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}_phyloFlash.report.csv",
        # phyloflash_fasta=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.all.final.fasta"

    conda:
        "../envs/phyloflash.yaml"

    params:
        db=config["phyloflash_db"],
        phylo_dir=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/"
    threads:
        config['threads']
    resources:
        mem_mb=config['max_mem']

    shell: # We have to zip the phyloflash output and move it because it will be stored in the WD otherwise.
        "phyloFlash.pl -dbhome {params.db} -lib {wildcards.sample_pool} -zip \
         -CPUs {threads} -read1 {input.forward_reads} -read2 {input.reverse_reads}; mv {wildcards.sample_pool}.phyloFlash.* {params.phylo_dir}"

rule unzip_phyloflash:
    input:
        phyloflash_tar=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.phyloFlash.tar.gz"
    output:
        phyloflash_output=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.all.final.fasta"
    params:
        phyloflash_dir=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/"
    shell:
        "tar -xzf {input.phyloflash_tar} -C {params.phyloflash_dir}"

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
        renamed_dir=f"{outdir}/results/07_maglinkage/{{sample_pool}}/"
    threads:
        config["threads"]
    resources:
        mem_mb=config['max_mem']
    shell:
        "MarkerMAG rename_reads -r1 {input.forward_reads} -r2 {input.reverse_reads} -p {wildcards.sample_pool} -fq \
        -t {threads}; mv {wildcards.sample_pool}_R*.fastq {params.renamed_dir}"

rule fastq_2_fasta:
    input:
        forward_reads=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R1.fastq",
        reverse_reads=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R2.fastq"
    output:
        fasta_forward=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R1.fasta",
        fasta_reverse=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R2.fasta"
    conda:
        "../envs/seqkit.yaml"
    threads:
        config["threads"]
    resources:
        mem_mb=config['max_mem']
    shell:
        """
        seqkit fq2fa {input.forward_reads} -o {output.fasta_forward} --threads {threads}
        seqkit fq2fa {input.reverse_reads} -o {output.fasta_reverse} --threads {threads}
        """


rule markermag_link:
    input:
        forward_reads=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R1.fasta",
        reverse_reads=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R2.fasta",
        phyloflash=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.all.final.fasta",
        mag_fasta=f"{outdir}/results/06_binning/drep/dereplicated_genomes",
        drep_done=f"{outdir}/results/06_binning/drep/dereplicated_genomes/drep.done"
    output:
        markerMAG_link=f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag/{{sample_pool}}_linkages_by_genome.txt",
        markerMAG_dir=directory(f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag"),
        markermag_done=touch(f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag/markermag.done")
    conda:
        "../envs/markerMAG.yaml"
    threads:
        config["threads"]
    resources:
        mem_mb=config['max_mem']
    shell:
        "MarkerMAG link -p {wildcards.sample_pool} -r1 {input.forward_reads} -r2 {input.reverse_reads} \
        -marker {input.phyloflash} -mag {input.mag_fasta} -o {output.markerMAG_dir} -x fa -t {threads} -force"
