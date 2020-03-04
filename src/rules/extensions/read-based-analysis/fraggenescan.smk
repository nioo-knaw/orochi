rule predict_genes:
    input:
        "scratch/filter/{sample}_R1.fasta"
    output:
        "scratch/read-based/fraggenescan/{sample}_forward_paired.ffn"
    params:
        prefix="scratch/read-based/fraggenescan/{sample}_forward_paired"
    conda:
        "../../envs/fraggenescan.yaml"
    threads: 32
    shell: "run_FragGeneScan.pl -genome={input} -out={params.prefix}  -complete=0  -train=illumina_5 -thread={threads}"
