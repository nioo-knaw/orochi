rule merge_and_rename:
    input:
        forward = lambda wildcards: config["data"][wildcards.sample]['forward'],
        reverse = lambda wildcards: config["data"][wildcards.sample]['reverse']
    output:
        forward=protected("scratch/unpack/{sample}_1.fastq.gz"),
        reverse=protected("scratch/unpack/{sample}_2.fastq.gz"),
    threads: 16
    run:
        if os.path.splitext(input[0])[1] == ".bz2":
            shell("pbzip2 -p{threads} -dc {input.forward} | pigz -p {threads} > {output.forward}")
            shell("pbzip2 -p{threads} -dc {input.reverse} | pigz -p {threads} > {output.reverse}")
        if os.path.splitext(input[0])[1] == ".gz":
            shell("pigz -p {threads} -dc {input.forward} | pigz -p {threads} > {output.forward}")
            shell("pigz -p {threads} -dc {input.reverse} | pigz -p {threads} > {output.reverse}")

