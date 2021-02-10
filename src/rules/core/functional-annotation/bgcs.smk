rule filter_contigs_antismash:
    input:
        "scratch/assembly/megahit/minimus2/primary.contigs.fa"
    output:
        "scratch/annotation/antismash/secondary.contigs.fa"
    params:
        length=5000
    conda:
        "../../../envs/seqtk.yaml"
    shell: 
       "seqtk seq -L {params.length} {input}  > {output}"

rule antismash:
    input:
        # TODO: Decide what is the input here
        #expand("scratch/assembly/megahit/{treatment}/{kmers}/final.contigs.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
        "scratch/annotation/antismash/secondary.contigs.fa"
    output:
        "results/annotation/antismash/secondary.contigs.gbk",
        "results/annotation/antismash/secondary.contigs.json"
    params:
        outdir="results/annotation/antismash/"
    conda:
        "../../../envs/antismash.yaml"
    threads: 16
    shell:
        "antismash --cb-general --cb-knownclusters --cb-subclusters --asf --genefinding-tool prodigal-m --output-dir {params.outdir} --cpus {threads} {input}"

rule test_antismash:
    input:
        "scratch/assembly/megahit/minimus2/primary.long.contigs.fa"
    output:
        "scratch/annotation/test_antismash/primary.long.contigs.gbk",
        "scratch/annotation/test_antismash/primary.long.contigs.json"
    params:
        outdir="scratch/annotation/test_antismash/"
    conda:
        "../../../envs/antismash.yaml"
    threads: 16
    shell:
        "antismash --cb-general --cb-knownclusters --cb-subclusters --asf --genefinding-tool prodigal-m --output-dir {params.outdir} --cpus {threads} {input}"

rule get_bgcs:
    input:
        "results/annotation/antismash/secondary.contigs.json"
    output:
        "scratch/annotation/antismash/bgcs.fasta"
    script:
        "../../../scripts/antismash_get_bgcs.py"

rule map_reads:
    input:
        bgcs="scratch/annotation/antismash/bgcs.fasta",
        forward=expand("scratch/host_filtering/{sample}_R1.fastq", project=config["project"], sample=config["data"]),
        reverse=expand("scratch/host_filtering/{sample}_R2.fastq", project=config["project"], sample=config["data"])
    output:
        "results/annotation/antismash/bgcs.count.txt"
    conda:
        "../../../envs/coverm.yaml"
    log: "scratch/annotation/antismash/bgcs.mapping.txt"
    threads: 24
    shell:
        "coverm contig --methods count --mapper minimap2-sr --proper-pairs-only -1 {input.forward} -2 {input.reverse} --reference {input.bgcs} --threads {threads} 2> {log} > {output}"

if config['big']=='bigscape':
   rule bigscape:
        input:
            gbks="results/annotation/antismash/secondary.contigs.gbk"
        params:
            inputdir="results/annotation/antismash"
        output:
            directory("results/annotation/bigscape")
        container:
            "docker://nselem/big-scape"
        conda:
            "../../../envs/bigscape.yaml"
        threads: 40
        shell:
            """
            git clone https://git.wur.nl/medema-group/BiG-SCAPE.git
            cd BiG-SCAPE
            wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz
            hmmpress Pfam-A.hmm
            python bigscape.py gbks -i {params.inputdir} -o {output} -c {threads}
            """

"""
if config['big']=='bigslice':
   rule bigslice:
        input: "scratch/annotation/antismash/secondary.contigs.gbk"
        params:
            inputdir="scratch/annotation/antismash"
            outdir="scratch/annotation/{big}"
        output:
        conda:
            "../../../envs/bigslice.yaml"
        log:
        threads:
        shell:
            "download_bigslice_hmmdb"
            "bigslice -i <{params.inputdir}> <{params.outdir}>"

rule rast:
# Uses myRAST batch processor to upload bigscape output (clusters) and gets rast IDs and downloads faa and txt files in batch.

rule prepare_corason:
# Arranges bigscape and rast output into specified input director structure needed for corason.

rule corason:
    input:
    output:
    container:
        "docker://nselem/corason"
    log:
    threads:
    shell:
"""
