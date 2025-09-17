# Rule to test the unfinished pipeline

rule downstream_test:
    input:
        # f"{outdir}/results/06_binning/metabat2/{{sample_pool}}/{{sample_pool}}_metabat2.done",
        # f"{outdir}/results/06_binning/maxbin2/{{sample_pool}}/{{sample_pool}}_maxbin2.done",
        # f"{outdir}/results/06_binning/dastool/{{sample_pool}}/{{sample_pool}}_DASTool.done",
        f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}/{{sample}}_assembly.fasta.gz",
        # f"{outdir}/results/06_binning/drep/dereplicated_genomes/drep.done",
        # test_target,
        # f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
        # f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta.gz",
        # f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_rev.fastq.gz",
        f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}/{{sample}}_orfs.fna",
        f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample}}/{{sample}}.contig2classification.names.summarise.txt",
        f"{outdir}/results/05_prokaryote_annotation/eggnog/{{sample}}/{{sample}}.emapper.annotations",
        # f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/temp_MetaPhlAn/{{sample_pool}}.txt",
        f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/merged_abundance_table.txt",
        f"{outdir}/results/03_assembly/size_filtered/{{sample}}_{minsize}/contigs_{{sample}}_{minsize}.fasta",
        f"{outdir}/results/04_gene_prediction/prodigal/{{sample}}/{{sample}}_proteins.faa",
        f"{outdir}/results/04_gene_prediction/whokaryote/{{sample}}/eukaryotes.fasta",
        f"{outdir}/results/04_gene_prediction/augustify/{{sample}}/{{sample}}_eukproteins.gff",
        f"{outdir}/results/06_binning/checkm2/{{sample}}/quality_report.tsv",
        f"{outdir}/results/06_binning/drep/checkm2_genomeinfo/{{sample}}_genomeinfo.tsv",
        # f"{outdir}/results/06_binning/drep/combined_genomeinfo.tsv",
        # f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag/{{sample_pool}}_linkages_by_genome.txt",
        # f"{outdir}/results/06_binning/BAT/{{sample_pool}}/{{sample_pool}}_BAT.done",
        # f"{outdir}/results/08_plots/{{sample_pool}}/{{sample_pool}}_krona.html",
        # f"{outdir}/results/08_plots/{{sample_pool}}/{{sample_pool}}_bins_scatterplot.png"
        # f"{outdir}/results/06_binning/BAT/{{sample_pool}}/{{sample_pool}}.bin2classification.names.txt"


    output:
#        test_file=f"{outdir}/results/05_test/{{sample}}/{{sample}}_test.txt",
        test_file1=f"{outdir}/results/05_test/{{sample}}/{{sample}}_test.txt"
    run:
        shell("echo {input}")
#        shell("touch {output.test_file} && touch {output.test_file1}")
        shell("touch {output.test_file1}")

        # shell("touch {output.test_file}")
