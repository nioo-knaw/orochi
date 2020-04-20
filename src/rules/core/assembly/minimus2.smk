rule merge_assemblies:
    input:
        expand("scratch/assembly/megahit/{treatment}/{kmers}/final.contigs.fa",treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        "scratch/assembly/megahit/primary.contigs.fa"
    shell:
        "cat {input} > {output}"

rule filter_contigs:
    input:
        "scratch/assembly/megahit/primary.contigs.fa"
    output:
        "scratch/assembly/megahit/primary.long.contigs.fa"
    params:
        length=100
    conda:
        "../../../envs/seqtk.yaml"
    shell: 
       "seqtk seq -L {params.length} {input}  > {output}"

rule contig_overlap:
    input:
        "scratch/assembly/megahit/primary.long.contigs.fa"
    output:
        "scratch/assembly/megahit/primary.long.contigs.99.fa"
    conda:
        "../../../envs/cd-hit.yaml"
    log:
       "scratch/assembly/megahit/cd-hit.log"
    shell:
        "cd-hit-est -i {input} -o {output} -T 90 -M 500000 -c 0.99 -n 10 > {log}"

rule contig_rename:
    input:
        "scratch/assembly/megahit/primary.long.contigs.fa"
    output:
        "scratch/assembly/megahit/primary.long.contigs.99.renamed.fa"
    shell:
        """awk '/^>/ {{print ">contig_" ++i; next}}{{print}}' < {input} > {output}"""

rule toAmos:
    input:
        "scratch/assembly/megahit/primary.long.contigs.99.renamed.fa"
    output:
        "scratch/assembly/megahit/primary.long.contigs.99.renamed.afg"
    conda:
        "../../../envs/amos.yaml"
    shell:
        "toAmos -s {input} -o {output}"

rule minimus2:
    input:
        "scratch/assembly/megahit/primary.long.contigs.99.renamed.afg"
    output:
        "scratch/assembly/megahit/primary.long.contigs.99.renamed.fasta",
        "scratch/assembly/megahit/primary.long.contigs.99.renamed.singletons.seq"
    conda:
        "../../../envs/amos.yaml"
    log:
        "scratch/assembly/megahit/primary.long.contigs.99.renamed.runAmos.log"
    shell:
        "minimus2 `file={input}; echo ${{file%.*}}` -D OVERLAP=100 MINID=95"

rule minimus2_merge:
    input:
        "scratch/assembly/megahit/primary.long.contigs.99.renamed.fasta",
        "scratch/assembly/megahit/primary.long.contigs.99.renamed.singletons.seq"
    output:  
        "scratch/assembly/megahit/secondary.contigs.fasta"
    shell:
        "cat {input} > {output}"
