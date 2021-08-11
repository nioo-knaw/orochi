rule merge_assemblies:
    input:
        expand("scratch/assembly/megahit/{treatment}/{kmers}/final.contigs.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        "scratch/assembly/megahit/minimus2/primary.contigs.fa"
    log:
        "logs/assembly/megahit/minimus2/merge_assemblies.log"
    conda:
        "../../../envs/orochi-base.yaml"
    shell:
        "cat {input} > {output} 2> {log}"

rule save_smalls:
    input: 
        "scratch/assembly/megahit/minimus2/primary.contigs.fa"
    output:
        "scratch/assembly/megahit/minimus2/primary.short.contigs.fa"
    log:
        "logs/assembly/megahit/minimus2/save_smalls.log"
    conda:
        "../../../envs/orochi-base.yaml"
    shell:
        """
        awk -v RS='>[^\\n]+\\n' 'length() <= 2000 {{printf "%s", prt $0}} {{prt = RT}}' {input} > {output} 2> {log}
        """

rule filter_contigs:
    input:
        "scratch/assembly/megahit/minimus2/primary.contigs.fa"
    output:
        "scratch/assembly/megahit/minimus2/primary.long.contigs.fa"
    params:
        length=config["filter_contigs_length"]
    log:
        "logs/assembly/megahit/minimus2/filter_contigs.log"
    conda:
        "../../../envs/seqtk.yaml"
    shell: 
       "seqtk seq -L {params.length} {input}  > {output} 2> {log}"

rule contig_overlap:
    input:
        "scratch/assembly/megahit/minimus2/primary.long.contigs.fa"
    output:
        "scratch/assembly/megahit/minimus2/primary.long.contigs.99.fa"
    conda:
        "../../../envs/cd-hit.yaml"
    log:
       "logs/assembly/megahit/minimus2/cd-hit.log"
    conda:
        "../../../envs/orochi-base.yaml"
    shell:
        "cd-hit-est -i {input} -o {output} -T 90 -M 500000 -c 0.99 -n 10 > {log}"

rule contig_rename:
    input:
        "scratch/assembly/megahit/minimus2/primary.long.contigs.fa"
    output:
        "scratch/assembly/megahit/minimus2/primary.long.contigs.99.renamed.fa"
    log:
        "logs/assembly/megahit/minimus2/contig_rename.log"
    conda:
        "../../../envs/orochi-base.yaml"
    shell:
        #"""awk '/^>/ {{print ">contig_" ++i; next}}{{print}}' < {input} > {output} 2> {log}"""
        "seqtk seq -C {input} | seqtk rename - contig_ > {output} 2> {log}"

rule toAmos:
    input:
        "scratch/assembly/megahit/minimus2/primary.long.contigs.99.renamed.fa"
    output:
        "scratch/assembly/megahit/minimus2/primary.long.contigs.99.renamed.afg"
    conda:
        "../../../envs/amos.yaml"
    log:
        "logs/assembly/megahit/minimus2/toAmos.log"
    shell:
        "toAmos -s {input} -o {output} 2> {log}"

rule minimus2:
    input:
        "scratch/assembly/megahit/minimus2/primary.long.contigs.99.renamed.afg"
    output:
        "scratch/assembly/megahit/minimus2/primary.long.contigs.99.renamed.fasta",
        "scratch/assembly/megahit/minimus2/primary.long.contigs.99.renamed.singletons.seq"
    conda:
        "../../../envs/amos.yaml"
    log:
        "logs/assembly/megahit/minimus2/primary.long.contigs.99.renamed.runAmos.log"
    shell:
        "minimus2 `file={input}; echo ${{file%.*}}` -D OVERLAP=100 MINID=95"

rule minimus2_merge:
    input:
        "scratch/assembly/megahit/minimus2/primary.long.contigs.99.renamed.fasta",
        "scratch/assembly/megahit/minimus2/primary.long.contigs.99.renamed.singletons.seq"
    output:  
        "scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    log:
        "logs/assembly/megahit/minimus2/minimus2_merge.log"
    conda:
        "../../../envs/orochi-base.yaml"
    shell:
        "cat {input} > {output} 2> {log}"

rule readd_smalls:
    input:
        short="scratch/assembly/megahit/minimus2/primary.short.contigs.fa",
        merged="scratch/assembly/megahit/minimus2/secondary.contigs.fasta"
    output: "scratch/assembly/megahit/minimus2/all.merged.contigs.fasta"
    log:
        "logs/assembly/megahit/minimus2/readd_smalls.log"
    conda:
        "../../../envs/orochi-base.yaml"
    shell: "cat {input.short} {input.merged} > {output} 2> {log}"
