rule das:
    input:
        #concoct="",
        #groopm="",
        #maxbin="",
        metabat="results/binning/metabat/metabat.bins.txt",
        vamb="results/binning/vamb/clusters.tsv",
        contigs="scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output: 
        "results/binning/DAS_tool/DASTool_summary.txt",
        "results/binning/DAS_tool/bin1.fasta"
    params:
        #outprefix="results/binning/DAS_tool/DAS"
        outprefix=lambda wildcards, output: output[0][:-16]
    log:
        "logs/binning/DAS_tool.log"
    conda:
        "../../../envs/dastool.yaml"
    shell:
        """
        DAS_Tool -i {input.metabat},
                    {input.vamb}
                 -l metabat,vamb
                 -c {input.contigs}
                 -o {params.outprefix}
                 --write_bins 1
                 2> {log}
        """
