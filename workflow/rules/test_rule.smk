# Rule to test the unfinished pipeline

rule downstream_test:
    input:
        # get_metabat_bins,
        # get_maxbin_bins,
        # get_dastool_bins,
        # get_drep_bins,
        # test_target,
        # 'config/samples_unsupervised.tsv',
        # rules.clustering.output,
        f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
        f"{outdir}/size_filter/contigs_{minsize}_{{sample_pool}}.fasta",
        f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_proteins.faa",
        f"{outdir}/results/04_gene_prediction/whokaryote/{{sample_pool}}/eukaryotes.fasta",
        f"{outdir}/results/04_gene_prediction/augustify/{{sample_pool}}/{{sample_pool}}_eukproteins.faa",

    output:
        test_file=f"{outdir}/results/05_test/{{sample_pool}}/{{sample_pool}}_test.txt"
    run:
        shell("echo {input}")
        shell("touch {output.test_file}")