# Rule to test the unfinished pipeline

rule downstream_test:
    input:
        # get_metabat_bins,
        # get_maxbin_bins,
        # get_dastool_bins,
        # get_drep_bins,
        # test_target,
        #'config/samples_unsupervised.tsv',
        #rules.clustering.output,
        f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",
        f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_orfs.fna",
        f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/{{sample_pool}}.contig2classification.names.summarise.txt",
        f"{outdir}/results/05_prokaryote_annotation/eggnog/{{sample_pool}}/{{sample_pool}}.emapper.annotations",
#        f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/{{sample}}.txt"
#        expand(f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/{samples}.txt", samples = SAMPLES)

    output:
#        test_file=f"{outdir}/results/05_test/{{sample}}/{{sample}}_test.txt",
        test_file1=f"{outdir}/results/05_test/{{sample_pool}}/{{sample_pool}}_test.txt"
    run:
        shell("echo {input}")
#        shell("touch {output.test_file} && touch {output.test_file1}")
        shell("touch {output.test_file1}")
