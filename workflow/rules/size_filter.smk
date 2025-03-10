""" Rules related to size filtering of contigs, gffs, and protein fasta files """

rule size_filter_contigs:
    input:
        contigs=branch(config['assembly_method'] == "coassembly",
            then=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample}}/{{sample}}_assembly.fasta",
            otherwise=f"{outdir}/results/03_assembly/single_sample_assembly/{{sample}}/{{sample}}_assembly.fasta")
    output:
        f"{outdir}/results/03_assembly/size_filtered/{{sample}}_{minsize}/contigs_{{sample}}_{minsize}.fasta"
    conda:
        "../envs/size_filter.yaml"
    params:
        size=config['min_contig_length'],
        outdir=f"{outdir}/results/03_assembly/size_filtered/{{sample}}_{minsize}",
        script=os.path.abspath("workflow/scripts/sizefilter_contigs.py")
    shell:
        "python {params.script} --contigs {input} --outdir {params.outdir} --minsize {params.size} --sample_name {wildcards.sample}"
