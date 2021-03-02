"""
rule vamb_filter:
    input: "scratch/assembly/megahit/{treatment}/{kmers}/final.contigs.fa"
    output: "scratch/vamb/contigs/{treatment}/{kmers}/long.contigs.fa"
    params:
        length=2000
    conda: "../../../envs/seqtk.yaml"
    shell: "seqtk seq -L {params.length} {input}  > {output}"
"""

rule concatenate:
    input: "scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output: "scratch/binning/vamb/catalogue.fna.gz"
    conda: "../../../envs/vamb.yaml"
    shell:     
        "concatenate.py {output} {input}"
    #Just cat with extras to make it more suitable to VAMB

rule vamb:
    input:
        catalogue="scratch/binning/vamb/catalogue.fna.gz"
        bam=expand("scratch/coverm/bamfiles/readsorted/{sample}.bam", sample=config["data"])
    output: 
        "results/binning/vamb/clusters.tsv",
        "results/binning/vamb/latent.npz",
        "results/binning/vamb/lengths.npz",
        "results/binning/vamb/log.txt",
        "results/binning/vamb/model.pt",
        "results/binning/vamb/mask.npz",
        "results/binning/vamb/tnf.npz"
    conda: "../../../envs/vamb.yaml"
    shell: 
        "rm -rf results/binning/vamb;"
        "vamb --outdir results/binning/vamb --fasta {input.catalogue} --bamfiles scratch/coverm/bamfiles/readsorted/*.bam -o C --minfasta 200000"

rule vamb_write_bins:
    input:
        clusters="results/binning/vamb/clusters.tsv",
        contigs="scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output: "results/binning/vamb/bins/bin1.fasta"
    params:
        outdir="results/binning/vamb"
    run:
        with open('{input.clusters}', 'w') as file:
            vamb.cluster.write_clusters(file, filtered_bins)
        keptcontigs = set.union(*filtered_bins.values())
        with open('{input.contigs}', 'rb') as file:
            fastadict = vamb.vambtools.loadfasta(file, keep=keptcontigs)
        bindir = '{params.outdir}'
        vamb.vambtools.write_bins(bindir, filtered_bins, fastadict, maxbins=500)
