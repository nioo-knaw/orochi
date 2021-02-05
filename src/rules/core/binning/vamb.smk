rule sample_assembly:
    input:         
        forward=lambda wildcards: expand("scratch/host_filtering/{sample}_R1.fastq", sample=config["data"]) if config['host_removal'] \
             else expand("scratch/filter/{sample}_R1.fasta", sample=config["data"]),
        reverse=lambda wildcards: expand("scratch/host_filtering/{sample}_R2.fastq", sample=config["data"]) if config['host_removal'] \
             else expand("scratch/filter/{sample}_R2.fasta", sample=config["data"]),
    output: temp("scratch/vamb/assembly/{sample}/{kmers}/contigs.fasta")
    params:
        outdir="scratch/vamb/assembly/{sample}/{kmers}",
        kmers = lambda wildcards: config["assembly-klist"][wildcards.kmers]
    log:
        "scratch/vamb/assembly/{sample}/{kmers}/spades.log"
    threads: 32
    conda: "../../../envs/spades.yaml"
    shell: "metaspades.py -m 1200 -1 {input.forward} -2 {input.reverse} --only-assembler -k {params.kmers} -t {threads} -o {params.outdir} --tmp-dir {params.outdir}/tmp/ 2>&1 > /dev/null"

rule vamb_filter_contigs:
    input:
        "scratch/vamb/assembly/{sample}/{kmers}/contigs.fasta"
    output:
        "scratch/vamb/assembly/{sample}/{kmers}/long.contigs.fasta"
    params:
        length=2000
    conda:
        "../../../envs/seqtk.yaml"
    shell: 
       "seqtk seq -L {params.length} {input}  > {output}"

rule concatenate:
    input: expand("scratch/vamb/assembly/{sample}/{kmers}/long.contigs.fasta", sample=config["data"])
    output: "scratch/vamb/catalogue.fna.gz"
    conda: "../../../envs/vamb.yaml"
    shell: "concatenate.py {output} {input}"

rule read_mapper:
    input:
        forward=lambda wildcards: expand("scratch/host_filtering/{sample}_R1.fastq", sample=config["data"]) if config['host_removal'] \
             else expand("scratch/filter/{sample}_R1.fasta", sample=config["data"]),
        reverse=lambda wildcards: expand("scratch/host_filtering/{sample}_R2.fastq", sample=config["data"]) if config['host_removal'] \
             else expand("scratch/filter/{sample}_R2.fasta", sample=config["data"]),
        catalogue="scratch/vamb/catalogue.fna.gz"
    output: "scratch/vamb/bamfiles/{sample}.bam"
    conda: "../../../envs/minimap2.yaml"
    shell:
        """
        minimap2 -d catalogue.mmi {input.catalogue}; # make index
minimap2 -t 8 -N 50 -ax sr catalogue.mmi {input.forward} {input.reverse} | samtools view -F 3584 -b --threads 8 > {output}
        """

rule vamb:
    input:
        catalogue="scratch/vamb/catalogue.fna.gz",
        bam=expand("scratch/vamb/bamfiles/{sample}.bam", sample=config["data"])
    output: "results/binning/vamb/clusters.tsv"
    params:
        outdir="results/binning/vamb"
    conda: "../../../envs/vamb.yaml"
    shell: "vamb --outdir {params.outdir} --fasta {input.catalogue} --bamfiles {input.bam} -o C --minfasta 200000"

"""
rule vamb_write_bins:
    input:
    output:
    conda: "../../../envs/vamb.yaml"
    shell:
"""
