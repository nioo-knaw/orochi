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
        catalogue="scratch/assembly/megahit/minimus2/secondary.contigs.fasta",
        bam=expand("scratch/coverm/bamfiles/readsorted/{sample}.bam", sample=config["data"])
    output: 
        "vamb/clusters.tsv",
        "vamb/latent.npz",
        "vamb/lengths.npz",
        "vamb/log.txt",
        "vamb/model.pt",
        "vamb/mask.npz",
        "vamb/tnf.npz"
    conda: "../../../envs/vamb.yaml"
    shell: "vamb --outdir vamb --fasta {input.catalogue} --bamfiles scratch/coverm/bamfiles/readsorted/*.bam -o C --minfasta 200000"

rule vamb_write_bins:
    input:
        clusters="vamb/clusters.tsv",
        contigs="scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output: "vamb/bins/bin1.fasta"
    params:
        outdir="vamb"
    run:
        with open('{input.clusters}', 'w') as file:
            vamb.cluster.write_clusters(file, filtered_bins)
        keptcontigs = set.union(*filtered_bins.values())
        with open('{input.contigs}', 'rb') as file:
            fastadict = vamb.vambtools.loadfasta(file, keep=keptcontigs)
        bindir = '{params.outdir}'
        vamb.vambtools.write_bins(bindir, filtered_bins, fastadict, maxbins=500)
