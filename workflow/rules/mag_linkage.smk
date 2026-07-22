""" Rules related to reconstructing 16S rRNA gene sequences and linking them to MAGs"""

phyloflash_version = config["phyloflash_version"]
phyloflash_dir = config["phyloflash_dir"]

# phyloFlash.pl locates its bundled PhyloFlash.pm/barrnap-HGV via FindBin, so it must be
# invoked from its extracted source tree rather than installed as a conda package.
rule download_phyloflash:
    output:
        script=f"{phyloflash_dir}/phyloFlash.pl"
    params:
        version=phyloflash_version,
        dir=phyloflash_dir
    log:
        f"{outdir}/logs/phyloflash/download_phyloflash.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.dir}
        curl -L https://github.com/HRGV/phyloFlash/archive/refs/tags/pf{params.version}.tar.gz \
            | tar -xz -C {params.dir} --strip-components=1 > {log} 2>&1
        find {params.dir} -type f -name "*.pl" -exec chmod +x {{}} + 2>> {log}
        find {params.dir} -type f -path "*/bin/*" -exec chmod +x {{}} + 2>> {log}
        find {params.dir} -type f -path "*/binaries/*" -exec chmod +x {{}} + 2>> {log}
        """

# We use the non-normalized reads (if coassembly) because the MAG coverage is also based on non-normalized reads
rule phyloflash:
    input:
        forward_reads = branch(config['assembly_method'] == "coassembly",
            then=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
            otherwise=f"{outdir}/results/02_filtered_reads/{{sample_pool}}_filt_1.fastq.gz"),
        reverse_reads = branch(config['assembly_method'] == "coassembly",
            then=f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_rev.fastq.gz",
            otherwise=f"{outdir}/results/02_filtered_reads/{{sample_pool}}_filt_2.fastq.gz"),
        phyloflash_script=f"{phyloflash_dir}/phyloFlash.pl"
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
    log:
        f"{outdir}/logs/phyloflash/phyloflash_{{sample_pool}}.log"

    shell: # We have to zip the phyloflash output and move it because it will be stored in the WD otherwise.
        """
        set -euo pipefail
        
        echo "Starting phyloFlash analysis for {wildcards.sample_pool}" > {log}
        echo "Timestamp: $(date)" >> {log}
        echo "Database: {params.db}" >> {log}
        echo "Forward reads: {input.forward_reads}" >> {log}
        echo "Reverse reads: {input.reverse_reads}" >> {log}
        echo "CPUs: {threads}" >> {log}
        echo "" >> {log}
        
        # Ensure output directory exists
        mkdir -p {params.phylo_dir}
        
        # Run phyloFlash
        {input.phyloflash_script} \
            -dbhome {params.db} \
            -lib {wildcards.sample_pool} \
            -zip \
            -CPUs {threads} \
            -read1 {input.forward_reads} \
            -read2 {input.reverse_reads} \
            >> {log} 2>&1
        
        # Move all phyloFlash output files to the target directory
        echo "" >> {log}
        echo "Moving phyloFlash output files to {params.phylo_dir}" >> {log}
        mv {wildcards.sample_pool}.phyloFlash.* {params.phylo_dir} 2>> {log}
        
        echo "phyloFlash analysis completed successfully at $(date)" >> {log}
        """

rule unzip_phyloflash:
    input:
        phyloflash_tar=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.phyloFlash.tar.gz"
    output:
        phyloflash_output=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.all.final.fasta",
        phyloflash_classification=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/{{sample_pool}}.phyloFlash.extractedSSUclassifications.csv"
    params:
        phyloflash_dir=f"{outdir}/results/07_maglinkage/{{sample_pool}}/phyloflash/"
    log:
        f"{outdir}/logs/phyloflash/unzip_phyloflash_{{sample_pool}}.log"
    shell:
        "tar -xzf {input.phyloflash_tar} -C {params.phyloflash_dir} 2> {log}"

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
    log:
        f"{outdir}/logs/unzip_reads/unzip_reads_{{sample_pool}}.log"
    shell:
        """
        set -euo pipefail
        
        echo "Starting read decompression for {wildcards.sample_pool}" > {log}
        echo "Timestamp: $(date)" >> {log}
        echo "Forward reads: {input.forward_reads}" >> {log}
        echo "Reverse reads: {input.reverse_reads}" >> {log}
        echo "" >> {log}
        
        # Decompress forward reads
        echo "Decompressing forward reads..." >> {log}
        gunzip -c {input.forward_reads} > {output.forward_out} 2>> {log}
        
        # Decompress reverse reads
        echo "Decompressing reverse reads..." >> {log}
        gunzip -c {input.reverse_reads} > {output.reverse_out} 2>> {log}
        
        # Verify output files
        forward_size=$(stat -f%z {output.forward_out} 2>/dev/null || stat -c%s {output.forward_out})
        reverse_size=$(stat -f%z {output.reverse_out} 2>/dev/null || stat -c%s {output.reverse_out})
        
        echo "" >> {log}
        echo "Decompression completed successfully at $(date)" >> {log}
        echo "Forward output size: $forward_size bytes" >> {log}
        echo "Reverse output size: $reverse_size bytes" >> {log}
        """

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
    log:
        f"{outdir}/logs/rename_reads/rename_reads_{{sample_pool}}.log"
    shell:
        """
        set -euo pipefail
        
        echo "Starting read renaming for {wildcards.sample_pool}" > {log}
        echo "Timestamp: $(date)" >> {log}
        echo "Input forward: {input.forward_reads}" >> {log}
        echo "Input reverse: {input.reverse_reads}" >> {log}
        echo "Threads: {threads}" >> {log}
        echo "" >> {log}
        
        # Run MarkerMAG rename_reads
        echo "Running MarkerMAG rename_reads..." >> {log}
        MarkerMAG rename_reads \
            -r1 {input.forward_reads} \
            -r2 {input.reverse_reads} \
            -p {wildcards.sample_pool} \
            -fq \
            -t {threads} \
            >> {log} 2>&1
        
        # Move renamed files to output directory
        echo "" >> {log}
        echo "Moving renamed files to {params.renamed_dir}" >> {log}
        mv {wildcards.sample_pool}_R*.fastq {params.renamed_dir} 2>> {log}
        
        # Verify output files exist
        if [ -f {output.forward_renamed} ] && [ -f {output.reverse_renamed} ]; then
            echo "Read renaming completed successfully at $(date)" >> {log}
            echo "Output files:" >> {log}
            echo "  - {output.forward_renamed}" >> {log}
            echo "  - {output.reverse_renamed}" >> {log}
        else
            echo "ERROR: Expected output files not found" >&2
            exit 1
        fi
        """

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
    log:
        f"{outdir}/logs/fastq_2_fasta/fastq_2_fasta_{{sample_pool}}.log"
    shell:
        """
        echo "Starting conversion of FASTQ to FASTA for {wildcards.sample_pool}" > {log}
        echo "Timestamp: $(date)" >> {log}
        echo "Input forward: {input.forward_reads}" >> {log}
        echo "Input reverse: {input.reverse_reads}" >> {log}
        echo "" >> {log}
        seqkit fq2fa {input.forward_reads} -o {output.fasta_forward} --threads {threads} 2>> {log}
        seqkit fq2fa {input.reverse_reads} -o {output.fasta_reverse} --threads {threads} 2>> {log}
        echo "Conversion completed successfully at $(date)" >> {log}
        """

rule build_markermag_mag_dir:
    input:
        drep_done=f"{outdir}/results/06_binning/drep/drep.done",
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
    log:
        f"{outdir}/logs/build_markermag_mag_dir/build_markermag_mag_dir_{{sample_pool}}.log"
    shell:
        r"""
        rm -rf {output.mag_dir}
        mkdir -p {output.mag_dir}

        python {params.script} \
            --drep-dir {params.drep_dir} \
            --input-bins {input.input_bins} \
            --sample-pool {wildcards.sample_pool} \
            --out-dir {output.mag_dir} \
            --copy-mode {params.copy_mode} \
            --done {output.done} \
            > {log} 2>&1
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
        markermag_done=touch(f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag/markermag.done")
    conda:
        "../envs/markerMAG.yaml"
    threads:
        min(config["threads"], 16)
    resources:
        mem_mb=config["max_mem"]
    log:
        f"{outdir}/logs/markermag_link/markermag_link_{{sample_pool}}.log"
    params:
        markermag_dir=f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag"
    shell:
        r"""
        set -euo pipefail
        
        # Initialize log
        echo "Starting MarkerMAG linkage analysis for {wildcards.sample_pool}" > {log}
        echo "Timestamp: $(date)" >> {log}
        echo "Threads: {threads}" >> {log}
        echo "" >> {log}
        
        # Log input files
        echo "Input files:" >> {log}
        echo "  Forward reads: {input.forward_reads}" >> {log}
        echo "  Reverse reads: {input.reverse_reads}" >> {log}
        echo "  PhyloFlash markers: {input.phyloflash}" >> {log}
        echo "  MAG directory: {input.mag_fasta}" >> {log}
        echo "" >> {log}
        
        # Validate and count MAGs
        echo "Validating MAG directory..." >> {log}
        n_mags=$(find -L {input.mag_fasta} -maxdepth 1 -type f \
            \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | wc -l)
        
        echo "Number of valid MAG files found: $n_mags" >> {log}
        
        if [ "$n_mags" -eq 0 ]; then
            echo "ERROR: No valid MAG fasta files found for {wildcards.sample_pool}" >> {log}
            echo "MAG directory contents:" >> {log}
            ls -lh {input.mag_fasta} >> {log} 2>&1 || echo "Directory not accessible" >> {log}
            exit 1
        fi
        
        # Check for broken symlinks
        echo "Checking for broken symlinks..." >> {log}
        broken_links=$(find {input.mag_fasta} -maxdepth 1 -xtype l 2>/dev/null | wc -l)
        
        if [ "$broken_links" -gt 0 ]; then
            echo "ERROR: Found $broken_links broken symlink(s) in {input.mag_fasta}" >> {log}
            echo "Broken symlinks:" >> {log}
            find {input.mag_fasta} -maxdepth 1 -xtype l -ls >> {log} 2>&1
            exit 1
        fi
        
        echo "MAG directory validation successful" >> {log}
        echo "" >> {log}
        
        # --- Dynamically size threads for MarkerMAG based on actual read count ---
        # Prevents handing MarkerMAG more threads than the data can meaningfully
        # use, which triggers an internal bug on very small datasets (empty
        # per-thread subsets during Rd2 -> missing merged file -> crash).
        n_reads=$(grep -c '^>' {input.forward_reads} || echo 0)
        echo "Read count in forward reads: $n_reads" >> {log}

        if [ "$n_reads" -lt 1000 ]; then
            effective_threads=1
        else
            effective_threads=$(( n_reads / 50 ))
            if [ "$effective_threads" -gt {threads} ]; then
                effective_threads={threads}
            fi
        fi
        echo "Effective threads passed to MarkerMAG: $effective_threads" >> {log}
        echo "" >> {log}
        
        # Prepare output directory
        echo "Preparing output directory: {params.markermag_dir}" >> {log}
        rm -rf {params.markermag_dir}
        mkdir -p {params.markermag_dir}
        
        # Run MarkerMAG link
        echo "Running MarkerMAG link analysis..." >> {log}
        echo "Command: MarkerMAG link -p {wildcards.sample_pool} -r1 {input.forward_reads} -r2 {input.reverse_reads} -marker {input.phyloflash} -mag {input.mag_fasta} -o {params.markermag_dir} -x fa -t {threads} -force" >> {log}
        echo "" >> {log}
        
        MarkerMAG link \
            -p {wildcards.sample_pool} \
            -r1 {input.forward_reads} \
            -r2 {input.reverse_reads} \
            -marker {input.phyloflash} \
            -mag {input.mag_fasta} \
            -o {params.markermag_dir} \
            -x fa \
            -t {threads} \
            -force \
            >> {log} 2>&1
        
        echo "" >> {log}
        echo "MarkerMAG link completed" >> {log}
        
        # Check if linkage file exists and has content
        if [ ! -f {output.markerMAG_link} ]; then
            echo "WARNING: MarkerMAG did not produce linkage output file" >> {log}
            echo "Creating empty placeholder with header" >> {log}
            printf "MarkerGene\tGenomicSeq\n" > {output.markerMAG_link}
        elif [ ! -s {output.markerMAG_link} ]; then
            echo "WARNING: MarkerMAG linkage file is empty" >> {log}
            echo "Adding header to empty file" >> {log}
            printf "MarkerGene\tGenomicSeq\n" > {output.markerMAG_link}
        else
            n_links=$(tail -n +2 {output.markerMAG_link} 2>/dev/null | wc -l)
            echo "Success: Found $n_links genome-marker linkages" >> {log}
        fi
        
        echo "" >> {log}
        echo "MarkerMAG linkage analysis completed successfully at $(date)" >> {log}
        
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
    log:
        f"{outdir}/logs/add_taxonomy_maglinkage/add_taxonomy_maglinkage_{{sample_pool}}.log"
    shell:
        """
        python3 workflow/scripts/add_markermag_taxonomy.py \
            -m {input.markermag_link} \
            -p {input.phyloflash_classification} \
            -o {output.tax_linked} \
            > {log} 2>&1
        """
