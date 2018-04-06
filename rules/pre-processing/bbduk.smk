rule filter_contaminants:
     input:
        forward="{project}/unpack/{sample}_1.fastq.gz",
        reverse="{project}/unpack/{sample}_2.fastq.gz",
     output:
        forward="{project}/filter/{sample}_R1.fasta",
        reverse="{project}/filter/{sample}_R2.fasta",
        stats="{project}/stats/{sample}_contaminants_stats.txt"
     params:
         phix="refs/phix.fasta",
         adapters="refs/illumina_scriptseq_and_truseq_adapters.fa",
         quality="25"
     log: "{project}/filter/{sample}.log"
     conda: "../../envs/bbmap.yaml"
     threads: 16
     shell:"""bbduk.sh -Xmx8g in={input.forward} in2={input.reverse} out={output.forward} out2={output.reverse} \
              ref={params.adapters},{params.phix} qtrim="rl" trimq={params.quality} threads={threads} stats={output.stats} 2> {log}"""

