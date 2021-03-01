"""
rule vamb_filter:
    input: "scratch/assembly/megahit/{treatment}/{kmers}/final.contigs.fa"
    output: "scratch/vamb/contigs/{treatment}/{kmers}/long.contigs.fa"
    params:
        length=2000
    conda: "../../../envs/seqtk.yaml"
    shell: "seqtk seq -L {params.length} {input}  > {output}"

rule concatenate:
    input: expand("scratch/vamb/contigs/{treatment}/{kmers}/long.contigs.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output: "results/binning/vamb/catalogue.fna.gz"
    conda: "../../../envs/vamb.yaml"
    shell:     
        "concatenate.py {output} {input}"
    #Just cat with extras to make it more suitable to VAMB
"""

rule vamb:
    input:
        catalogue="scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
        #bam="scratch/coverm/bamfiles/secondary.contigs.fasta.{sample}_R1.fastq.bam"
    output: "results/binning/vamb/clusters.tsv"
    conda: "../../../envs/vamb.yaml"
    shell: "vamb --outdir results/binning/vamb --fasta {input.catalogue} --bamfiles scratch/coverm/bamfiles/*.bam -o C --minfasta 200000"

rule vamb_write_bins:
    input:
        clusters="results/binning/vamb/clusters.tsv",
        contigs=expand("scratch/vamb/contigs/{treatment}/{kmers}/long.contigs.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output: "results/vamb/bins/bin1.fasta"
    params:
        outdir="results/vamb/bins"
    run:
        with open('{input.clusters}', 'w') as file:
            vamb.cluster.write_clusters(file, filtered_bins)
        keptcontigs = set.union(*filtered_bins.values())
        with open('{input.contigs}', 'rb') as file:
            fastadict = vamb.vambtools.loadfasta(file, keep=keptcontigs)
        bindir = '{params.outdir}'
        vamb.vambtools.write_bins(bindir, filtered_bins, fastadict, maxbins=500)
