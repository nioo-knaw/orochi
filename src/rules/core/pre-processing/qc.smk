rule filter:
     input:
        forward="scratch/unpack/{sample}_1.fastq",
        rev="scratch/unpack/{sample}_2.fastq"
     output:
        forward=temp("scratch/filter/{sample}_R1.fq"),
        rev=temp("scratch/filter/{sample}_R2.fq"),
        stats="scratch/stats/{sample}_contaminants_stats.txt"
     params:
         phix="refs/phix.fasta",
         adapters="refs/illumina_scriptseq_and_truseq_adapters.fa",
         ref = "src/refs/nextera.fa.gz",
         quality="25"
     log: "logs/filter/{sample}.log"
     conda: "../../../envs/bbmap.yaml"
     threads: 16
     shell:"""bbduk.sh in={input.forward} in2={input.rev} out={output.forward} out2={output.rev} \
     trimpolygright=1 \
     entropy=0.6 entropywindow=50 entropymask=f \
     qtrim=rl trimq={params.quality} \
     minlength=51 \
     ref={params.ref} ktrim=r \
     stats={output.stats} \
     t={threads} 2> {log}"""

rule phix_removal:
    input:
        forward="scratch/filter/{sample}_R1.fq",
        rev="scratch/filter/{sample}_R2.fq",
    output:
        forward=temp("scratch/filter/{sample}_R1.nophix.fq"),
        rev=temp("scratch/filter/{sample}_R2.nophix.fq"),
    params:
        ref = "src/refs/phix174_ill.ref.fa.gz",
    log: "logs/filter/phix_removal_{sample}.log"
    conda: "../../../envs/bbmap.yaml"
    threads: 16
    shell: "bbmap.sh ref={params.ref} in1={input.forward} in2={input.rev} outu1={output.forward} outu2={output.rev} t={threads} 2> {log}"

rule index_host:
    input:
        fasta=config["reference"]
    output:
        index=config["reference"] + ".bwt"
    log: "logs/filter/index_host.log"
    conda: "../../../envs/bwa.yaml"
    shell: "bwa index {input}"

rule map_to_host:
    input:
        forward="scratch/filter/{sample}_R1.nophix.fq",
        rev="scratch/filter/{sample}_R2.nophix.fq",
        index=config["reference"] + ".bwt"
    output:
        "scratch/host_filtering/{sample}.sam"
    params:
        refindex=lambda wildcards, input: os.path.splitext(input[2]) [0]
    conda: "../../../envs/bwa.yaml"
    log: "logs/filter/bwa_{sample}.log"
    threads: 16
    shell: "bwa mem -t {threads} {params.refindex} {input.forward} {input.rev} -o {output} 2> {log}"

rule mapping_stats:
    input:
        "scratch/host_filtering/{sample}.sam"
    output:
        "scratch/host_filtering/{sample}.flagstat.txt"
    conda: "../../../envs/samtools.yaml"
    log: "logs/filter/mapping_stats_{sample}.log"
    shell: "samtools flagstat {input} > {output} 2> {log}"

rule get_unmapped:
    input:
        "scratch/host_filtering/{sample}.sam"
    output:
        "scratch/host_filtering/{sample}.unmapped.bam"
    conda: "../../../envs/samtools.yaml"
    log : "logs/filter/get_unmapped_{sample}.log"
    # TODO: what is a good quality value? With -f 4 also alignments with q=0 are reported. For bwa this should mean mapping to multiple locations. From q=50 onwards a blastn also hits the reference genome.
    shell: "samtools view -b -f 4 {input} > {output} 2> {log}"

rule sort_unmapped:
    input:
        "scratch/host_filtering/{sample}.unmapped.bam"
    output:
        "scratch/host_filtering/{sample}.unmapped.sorted.bam"
    log: "logs/filter/sort_unmapped_{sample}.log"
    conda: "../../../envs/samtools.yaml"
    shell: "samtools sort {input} > {output}"

rule bamToFastq_unmapped:
    input:
        "scratch/host_filtering/{sample}.unmapped.sorted.bam"
    output:
        forward="scratch/host_filtering/{sample}_R1.fastq",
        rev="scratch/host_filtering/{sample}_R2.fastq"
    log: "logs/host_filtering/{sample}_bamtofastq.log"
    conda: "../../../envs/bedtools.yaml"
    shell: "bamToFastq -i {input}  -fq {output.forward} -fq2 {output.rev} 2> {log}"

rule get_mapped:
    input:
        "scratch/host_filtering/{sample}.sam"
    output:
        "scratch/host_filtering/{sample}.mapped.bam"
    conda: "../../../envs/samtools.yaml"
    log: "logs/filter/get_mapped_{sample}.log"
    shell: "samtools view -b -F 4 {input} > {output} 2> {log}"

rule sort_mapped:
    input:
        "scratch/host_filtering/{sample}.mapped.bam"
    output:
        "scratch/host_filtering/{sample}.mapped.sorted.bam"
    conda: "../../../envs/samtools.yaml"
    log: "logs/filter/sort_mapped_{sample}.log"
    shell: "samtools sort {input} > {output}"

rule bamToFastq_mapped:
    input:
        "scratch/host_filtering/{sample}.mapped.sorted.bam"
    output:
        forward="scratch/host_filtering/{sample}_R1.mapped.fastq",
        rev="scratch/host_filtering/{sample}_R2.mapped.fastq"
    log: "logs/host_filtering/{sample}_bamtofastq.mapped.log"
    conda: "../../../envs/bedtools.yaml"
    shell: "bamToFastq -i {input}  -fq {output.forward} -fq2 {output.rev} 2> {log}"

