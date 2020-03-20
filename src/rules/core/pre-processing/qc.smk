rule filter:
     input:
        forward="scratch/unpack/{sample}_1.fastq",
        reverse="scratch/unpack/{sample}_2.fastq"
     output:
        forward="scratch/filter/{sample}_R1.fasta",
        reverse="scratch/filter/{sample}_R2.fasta",
        stats="scratch/stats/{sample}_contaminants_stats.txt"
     params:
         phix="refs/phix.fasta",
         adapters="refs/illumina_scriptseq_and_truseq_adapters.fa",
         quality="25"
     log: "scratch/filter/{sample}.log"
     conda: "../../../envs/bbmap.yaml"
     threads: 16
     shell:"""bbduk.sh in={input.forward} in2={input.reverse} out={output.forward} out2={output.reverse} \
     trimpolygright=1 \
     entropy=0.6 entropywindow=50 entropymask=f \
     qtrim=rl trimq={params.quality} \
     minlength=51 \
     ref=$CONDA_PREFIX/opt/bbmap-38.79-0/resources/nextera.fa.gz ktrim=r \
     stats={output.stats} \
     t={threads} 2> {log}"""

rule host_removal:
    input:
        forward="scratch/filter/{sample}_R1.fasta",
        reverse="scratch/filter/{sample}_R2.fasta",
    output:
        forward="scratch/host_filtering/{sample}_R1.fasta",
        reverse="scratch/host_filtering/{sample}_R2.fasta",
    params:
        refindex=config["reference_index"],
        reference=config["reference"]
    log: "scratch/filter/host_removal_{sample}.log"
    conda: "../../../envs/bbmap.yaml"
    threads: 16
    shell: "bbsplit.sh ref=$CONDA_PREFIX/opt/bbmap-38.79-0/resources/phix174_ill.ref.fa.gz,{params.reference} in1={input.forward} in2={input.reverse} basename=%.fasta outu1={output.forward} outu2={output.reverse} t={threads}"



