rule prodigal:
    input:
        contigs=branch(config['assembly_method'] == "coassembly",
            then=f"{outdir}/results/03_assembly/coassembly/assembly_{{sample_pool}}/{{sample_pool}}_assembly.fasta",
            otherwise=f"{outdir}/results/03_assembly/single_sample_assembly/{{sample_pool}}/{{sample_pool}}_assembly.fasta")

    output:
        gff=f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_genes.gff",
        faa=f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_proteins.faa",
        fna=f"{outdir}/results/04_gene_prediction/prodigal/{{sample_pool}}/{{sample_pool}}_orfs.fna"

    conda:
        "../envs/prodigal.yaml"

    shell:
        "prodigal -i {input.contigs} -o {output.gff} -a {output.faa} -d {output.fna} -p meta -f gff"

rule dram:
    input: rules.prodigal.output.fna
    output: "{outdir}/results/05_prokaryotic_annotation"
    params:

    conda:
        "../envs/dram.yaml"
    shell:
        "DRAM.py annotate -i 'my_bins/*.fa' -o {output}"
