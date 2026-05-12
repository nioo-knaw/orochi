""" Rules related to reconstructing 16S rRNA gene sequences and linking them to MAGs"""

# We use the non-normalized reads (if coassembly) because the MAG coverage is also based on non-normalized reads
rule phyloflash:
    input:
        forward_reads = branch(config['assembly_method'] == "coassembly",
            then=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
            otherwise=f"{outdir}/results/02_filtered_reads/{{sample_pool}}_filt_1.fastq.gz"),
        reverse_reads = branch(config['assembly_method'] == "coassembly",
            then=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_rev.fastq.gz",
            otherwise=f"{outdir}/results/02_filtered_reads/{{sample_pool}}_filt_2.fastq.gz")
        # forward_reads=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
        # reverse_reads=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_rev.fastq.gz"

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
        phyloflash_output=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.all.final.fasta",
        phyloflash_classification=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.phyloFlash.extractedSSUclassifications.csv"
    params:
        phyloflash_dir=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/"
    shell:
        "tar -xzf {input.phyloflash_tar} -C {params.phyloflash_dir}"

rule unzip_reads:
    input:
        forward_reads = branch(config['assembly_method'] == "coassembly",
            then=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
            otherwise=f"{outdir}/results/02_filtered_reads/{{sample_pool}}_filt_1.fastq.gz"),
        reverse_reads = branch(config['assembly_method'] == "coassembly",
            then=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_rev.fastq.gz",
            otherwise=f"{outdir}/results/02_filtered_reads/{{sample_pool}}_filt_2.fastq.gz")
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

rule build_markermag_mag_dir:
    input:
        drep_done=f"{outdir}/results/06_binning/drep/dereplicated_genomes/drep.done",
        input_bins=f"{outdir}/results/06_binning/drep/input_bins.txt"
    output:
        mag_dir=directory(
            f"{outdir}/results/06_binning/drep/markermag_by_sample/{{sample_pool}}/dereplicated_genomes"
        ),
        done=touch(
            f"{outdir}/results/06_binning/drep/markermag_by_sample/{{sample_pool}}/dereplicated_genomes/.done"
        )
    conda:
        "../envs/drep.yaml"
    params:
        drep_dir=f"{outdir}/results/06_binning/drep",
        script=os.path.abspath("workflow/scripts/build_markermag_mag_dirs.py"),
        copy_mode=config.get("markermag_mag_copy_mode", "symlink")
    shell:
        """
        python {params.script} \
            --drep-dir {params.drep_dir} \
            --input-bins {input.input_bins} \
            --sample-pool {wildcards.sample_pool} \
            --out-dir {output.mag_dir} \
            --copy-mode {params.copy_mode} \
            --done {output.done}
        """

rule markermag_link:
    input:
        forward_reads=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R1.fasta",
        reverse_reads=f"{outdir}/results/07_maglinkage/{{sample_pool}}/{{sample_pool}}_R2.fasta",
        phyloflash=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.all.final.fasta",
        mag_fasta=f"{outdir}/results/06_binning/drep/markermag_by_sample/{{sample_pool}}/dereplicated_genomes",
        mag_fasta_done=f"{outdir}/results/06_binning/drep/markermag_by_sample/{{sample_pool}}/dereplicated_genomes/.done"
    output:
        markerMAG_link=f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag/{{sample_pool}}_linkages_by_genome.txt",
        markerMAG_dir=directory(f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag"),
        markermag_done=touch(f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag/markermag.done")
    conda:
        "../envs/markerMAG.yaml"
    threads:
        64
        # This is a maximum number of threads, MarkerMAG will use less if not available.
        # With more threads the file number limit will be reached (>1024 open files).
    resources:
        mem_mb=config['max_mem']
    log:
        f"{outdir}/logs/markermag/{{sample_pool}}.log"
    shell:
        r"""
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        set -x

        echo "[MarkerMAG] sample_pool={wildcards.sample_pool}"
        echo "[MarkerMAG] MAG directory: {input.mag_fasta}"
        echo "[MarkerMAG] Number of MAGs:"
        find {input.mag_fasta} -maxdepth 1 -type f -name "*.fa" | wc -l

        MarkerMAG link \
            -p {wildcards.sample_pool} \
            -r1 {input.forward_reads} \
            -r2 {input.reverse_reads} \
            -marker {input.phyloflash} \
            -mag {input.mag_fasta} \
            -o {output.markerMAG_dir} \
            -x fa \
            -t {threads} \
            -force

        touch {output.markermag_done}
        """

rule add_taxonomy_maglinkage:
    input:
        markermag_link=f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag/{{sample_pool}}_linkages_by_genome.txt",
        phyloflash_classification=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.phyloFlash.extractedSSUclassifications.csv"
    output:
        tax_linked=f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag/{{sample_pool}}_linkages_by_genome_taxonomy.txt"
    conda:
        "../envs/python_simple.yaml"
    threads:
        config['threads']
    resources:
        mem_mb=config['max_mem']
    shell:
        """
        python3 workflow/scripts/add_markermag_taxonomy.py -m {input.markermag_link} -p {input.phyloflash_classification} -o {output.tax_linked}
        """

