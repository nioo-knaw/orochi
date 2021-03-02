"""
rule concoct_cutup:
    input: "scratch/assembly/megahit/minimus2/secondary.contigs.fa"
    output: "scratch/binning/concoct/contigs_10K.fa"
    conda: "../../../envs/concoct.yaml"
    shell: "cut_up_fasta.py {input} -c 10000 -o 0 --merge_last -b contigs_10K.bed > {output}"

rule concoct_coverage:
    input: "scratch/assembly/megahit/minimus2/secondary.contigs.fa"
    output: "scratch/binning/concoct/contigs_10K.fa"
    conda: "../../../envs/concoct.yaml"
    shell: "cut_up_fasta.py {input} -c 10000 -o 0 --merge_last -b contigs_10K.bed > {output}"
"""
rule concoct:
    input:
        contigs="scratch/assembly/megahit/minimus2/secondary.contigs.fa",
        coverage="results/stats/coverage.tsv"
    output: "results/binning/concoct/clustering_merged.csv"
    params:
        outdir: "results/binning/concoct/"
    conda: "../../../envs/concoct.yaml"
    shell: "concoct --composition_file {input.contigs} --coverage_file {input.coverage} -b {params.outdir}"

rule concoct_write_bins:
    input:
        contigs="scratch/assembly/megahit/minimus2/secondary.contigs.fa",
        clusters="results/binning/concoct/clustering_merged.csv"
    output: "results/binning/concoct/fasta_bins/bin1.fasta"
    params:
        outdir: "results/binning/concoct/fasta_bins"
    conda: "../../../envs/concoct.yaml"
    shell:
        """
        extract_fasta_bins.py {input.contigs} {input.clusters} --output_path {params.outdir}
        """
