rule das:
    input:
        vamb="results/binning/vamb/clusters.tsv",
        metabat="results/binning/mmgenome/metabat.bins.txt",
        contigs="scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output: 
        "results/binning/DAS_tool/DASTool_summary.txt",
        "results/binning/DAS_tool/bin1.fasta"
    params:
        outprefix="results/binning/DAS_tool/DAS"
    conda:
        "../../../envs/dastool.yaml"
    shell:
        """
        DAS_tool -i {input.vamb},
                    {input.metabat}
                 -l vamb,metabat
                 -c {input.contigs}
                 -o {params.outprefix}
                 --write_bins 1
        """
