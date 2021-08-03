rule merge_and_rename:
    input:
        forward = lambda wildcards: config["data"][wildcards.sample]['forward'],
        rev = lambda wildcards: config["data"][wildcards.sample]['rev']
    output:
        forward=protected("scratch/unpack/{sample}_1.fastq"),
        rev=protected("scratch/unpack/{sample}_2.fastq"),
    conda: "../../../envs/pigz.yaml"
    log: 
        forward = "logs/data.input_{sample}_1.log",
        rev = "logs/data.input_{sample}_2.log"
    threads: 16
    shell:
        """
        input=({input.forward})
        if [ "${{input##*.}}" == "bz2" ]; then
            pbzip2 -p{threads} -dc {input.forward}  > {output.forward} 2> {log.forward}
            pbzip2 -p{threads} -dc {input.rev}  > {output.rev} 2> {log.rev}
        elif [ "${{input##*.}}" == "gz" ]; then
            pigz -p {threads} -dc {input.forward}  > {output.forward} 2> {log.forward}
            pigz -p {threads} -dc {input.rev}  > {output.rev} 2> {log.rev}
        else
            cp {input.forward} {output.forward}
            cp {input.rev} {output.rev}
        fi
        """
