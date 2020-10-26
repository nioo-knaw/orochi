rule merge_and_rename:
    input:
        forward = lambda wildcards: config["data"][wildcards.sample]['forward'],
        reverse = lambda wildcards: config["data"][wildcards.sample]['reverse']
    output:
        forward=protected("scratch/unpack/{sample}_1.fastq"),
        reverse=protected("scratch/unpack/{sample}_2.fastq"),
    conda: "../../../envs/pigz.yaml"
    threads: 16
    shell:
        if [{input.forward} = "*.bz2"]; then
            "pbzip2 -p{threads} -dc {input.forward}  > {output.forward}"
            "pbzip2 -p{threads} -dc {input.reverse}  > {output.reverse}"
        if [{input.forward} = "*.gz"]; then
            "pigz -p {threads} -dc {input.forward}  > {output.forward}"
            "pigz -p {threads} -dc {input.reverse}  > {output.reverse}"
        else
            "cp {input.forward} {output.forward}"
            "cp {input.reverse} {output.reverse}"

    """
    run:
        if os.path.splitext(input[0])[1] == ".bz2":
            shell("pbzip2 -p{threads} -dc {input.forward}  > {output.forward}")
            shell("pbzip2 -p{threads} -dc {input.reverse}  > {output.reverse}")
        if os.path.splitext(input[0])[1] == ".gz":
            shell("pigz -p {threads} -dc {input.forward}  > {output.forward}")
            shell("pigz -p {threads} -dc {input.reverse}  > {output.reverse}")
        else:
            shell("cp {input.forward} {output.forward}")
            shell("cp {input.reverse} {output.reverse}")
    """