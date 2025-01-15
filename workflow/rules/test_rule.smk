# Rule to test the unfinished pipeline

if config['sample_pooling'] == "supervised":
    test_target = expand(f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta.gz",
        sample_pool=samples["sample_pool"].unique())
elif config['sample_pooling'] == "unsupervised":
    test_target = rules.clustering.output
else:
    raise ValueError("Unsupported sample pooling option, please select between 'supervised' and 'unsupervised'")

rule downstream_test:
    input:
        # get_metabat_bins,
        # get_maxbin_bins,
        # get_dastool_bins,
        # get_drep_bins,
        # test_target,
        'config/samples_unsupervised.tsv',
        f"{outdir}/results/03_assembly/coassembly/pools/{{sample_pool}}_forward.fastq.gz",

    output:
        test_file=f"{outdir}/results/05_test/{{sample_pool}}/{{sample_pool}}_test.txt"
    run:
        shell("echo {input}")
        shell("touch {output.test_file}")