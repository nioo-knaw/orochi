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
    input: "scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output: "scratch/vamb/catalogue.fna.gz"
    conda: "../../../envs/vamb.yaml"
    shell: "concatenate.py {output} {input}"

rule vamb:
    input:
        catalogue="scratch/vamb/catalogue.fna.gz",
        bam=lambda wildcards: "scratch/coverm/bamfiles/megahit/{treatment}/{kmers}/assembly.fa.{sample}_R1.fastq.bam", sample=config["data"], treatment=config["treatment"], kmers=config["kmers"]
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
