rule predict_genes:
    input:
        "{project}/filter/{sample}_R1.fasta"
    output:
        "{project}/read-based/fraggenescan/{sample}_forward_paired"
    conda:
        "../../envs/fraggenescan.yaml"
    threads: 32
    shell: "run_FragGeneScan.pl -genome={input} -out={output}  -complete=0  -train=illumina_5 -thread={threads}"
