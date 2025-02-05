rule prodigal:
    input:
        contigs=branch(config['assembly_method'] == "coassembly",
            then=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta",
            otherwise=f"{outdir}/results/03_assembly/single_sample_assembly/{{sample_pool}}/{{sample_pool}}_assembly.fasta")
    output:
        gff = f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_genes.gff",
        faa = f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_proteins.faa",
        fna = f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_orfs.fna"
    conda:
        "../envs/prodigal.yaml"
    shell:
        "prodigal -i {input.contigs} -o {output.gff} -a {output.faa} -d {output.fna} -p meta -f gff"


rule CAT:
    input:
        contigs = {rules.megahit.output.contigs},
        proteins = {rules.prodigal.output.faa}
    output:
        CAT = f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/{{sample_pool}}.contig2classification.txt",
        names = f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/{{sample_pool}}.contig2classification.names.txt",
        summary = f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/{{sample_pool}}.contig2classification.names.summarise.txt"
    params:
        db = config["CAT_database"],
        tax = config["CAT_taxonomy"],
        out_dir = f"{outdir}/results/05_prokaryote_annotation/CAT/{{sample_pool}}/",
        prefix = f"{{sample_pool}}",
        threads = config["threads"]
    conda:
        "../envs/cat.yaml"
    shell:
        """
        mkdir -p {params.out_dir}
        CAT_pack contigs -c {input.contigs} -n {params.threads} -d {params.db} -t {params.tax} -p {input.proteins} -o {params.out_dir}{params.prefix}
        CAT_pack add_names -i {output.CAT} -o {output.names} -t {params.tax} --only_official --exclude_scores
        CAT_pack summarise -c {input.contigs} -i {output.names} -o {output.summary}
        """






#rule dram:
#    input: rules.prodigal.output.fna
#    output: "{outdir}/results/05_prokaryotic_annotation"
#    params:
#
#    conda:
#        "../envs/dram.yaml"
#    shell:
#        "DRAM.py annotate -i 'my_bins/*.fa' -o {output}"

