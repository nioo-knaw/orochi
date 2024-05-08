if config['assembly_method']=='single assembly':
    rule spades:
        input:
            forward = "results/02_filtered_reads/{sample}_filt_1.fastq.gz",
            rev = "results/02_filtered_reads/{sample}_filt_2.fastq.gz",
        output:
            "results/03_single_sample_assembly/{sample}/contigs.fasta",
        params:
            outdir = "results/03_single_sample_assembly/{sample}",
            kmers = config['kmers'],
        conda:
            "../envs/single_assembly.yaml"
        shell: "spades.py --meta -m 1200 -1 {input.forward} -2 {input.rev} --only-assembler -k {params.kmers} -t {threads} -o {params.outdir} --tmp-dir {params.outdir}/tmp/" 

    rule rename_spades:
        input:
            contigs = rules.spades.output
        output:
            gzip = protected("results/03_single_sample_assembly/{sample}/assembly.fasta.gz"),
            fasta = temp("results/03_single_sample_assembly/{sample}/assembly.fasta")
        run:
            shell("cat {input.contigs} | awk '{{print $1}}' | sed 's/NODE/contig/' > {output.fasta}")
            shell("gzip -c {output.fasta} > {output.gzip}")
