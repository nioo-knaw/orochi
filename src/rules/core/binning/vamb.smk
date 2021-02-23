rule vamb_filter:
    input: "scratch/assembly/megahit/{treatment}/{kmers}/final.contigs.fa"
    output: "scratch/vamb/contigs/{treatment}/{kmers}/long.contigs.fa"
    params:
        length=2000
    conda: "../../../envs/seqtk.yaml"
    shell: "seqtk seq -L {params.length} {input}  > {output}"

rule concatenate:
    input: expand("scratch/vamb/contigs/{treatment}/{kmers}/long.contigs.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output: "scratch/vamb/catalogue.fna.gz"
    conda: "../../../envs/vamb.yaml"
    shell: "concatenate.py {output} {input}"
    #Just cat with extras to make it more suitable to VAMB

rule vamb:
    input:
        catalogue="scratch/vamb/catalogue.fna.gz",
        bam=expand("scratch/coverm/bamfiles/{{treatment}}/{{kmers}}/assembly.fa.{sample}_R1.fastq.bam", sample=config["data"])
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
