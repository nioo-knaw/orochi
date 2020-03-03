rule filter_contaminants:
     input:
        forward=rules.merge_and_rename.output.forward,
        reverse=rules.merge_and_rename.output.reverse
     output:
        forward="scratch/filter/{sample}_R1.fasta",
        reverse="scratch/filter/{sample}_R2.fasta",
        stats="scratch/stats/{sample}_contaminants_stats.txt"
     params:
         phix="refs/phix.fasta",
         adapters="refs/illumina_scriptseq_and_truseq_adapters.fa",
         quality="25"
     log: "scratch/filter/{sample}.log"
     conda: "../../envs/bbmap.yaml"
     threads: 16
     shell:"""bbduk.sh -Xmx8g in={input.forward} in2={input.reverse} out={output.forward} out2={output.reverse} \
              ref={params.adapters},{params.phix} qtrim="rl" trimq={params.quality} threads={threads} stats={output.stats} 2> {log}"""

