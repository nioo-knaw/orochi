# Rule to test the unfinished pipeline

rule downstream_test:
    input:
        get_metabat_bins,
        get_maxbin_bins,
        get_dastool_bins,
        get_drep_bins,
        # test_target,
        f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
        f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta.gz",
        # f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_rev.fastq.gz",
        f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_orfs.fna",
        f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/{{sample_pool}}.contig2classification.names.summarise.txt",
        # f"{outdir}/results/05_prokaryote_annotation/eggnog/{{sample_pool}}/{{sample_pool}}.emapper.annotations",
        # f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/merged_abundance_table.tsv",
#        f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/{{sample}}.txt"
#        expand(f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/{samples}.txt", samples = SAMPLES)
        f"{outdir}/results/03_assembly/size_filtered/{{sample_pool}}_{minsize}/contigs_{{sample_pool}}_{minsize}.fasta",
        f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_proteins.faa",
        f"{outdir}/results/04_gene_prediction/whokaryote/{{sample_pool}}/eukaryotes.fasta",
        f"{outdir}/results/04_gene_prediction/augustify/{{sample_pool}}/{{sample_pool}}_eukproteins.gff",
        f"{outdir}/results/06_binning/checkm2/quality_report.tsv"


    output:
#        test_file=f"{outdir}/results/05_test/{{sample}}/{{sample}}_test.txt",
        test_file1=f"{outdir}/results/05_test/{{sample_pool}}/{{sample_pool}}_test.txt"
    run:
        shell("echo {input}")
#        shell("touch {output.test_file} && touch {output.test_file1}")
        shell("touch {output.test_file1}")

        # shell("touch {output.test_file}")
