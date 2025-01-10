# Rule to test the unfinished pipeline

rule downstream_test:
    input:
        # get_metabat_bins,
        # get_maxbin_bins,
        # get_dastool_bins,
        get_drep_bins

    output:
        test_file=f"{outdir}/results/05_test/{{sample_pool}}/{{sample_pool}}_test.txt"
    run:
        shell("echo {input}")
        shell("touch {output.test_file}")