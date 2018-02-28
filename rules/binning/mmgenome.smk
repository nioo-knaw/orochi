rule mmgenome_coverage:
    input:
        expand("{{project}}/bamm/{{assembler}}/{{treatment}}/{{kmers}}/assembly.{sample}_forward_paired.bam", sample=config["data"]) 
    output:
        "{project}/mmgenome/{assembler}/{treatment}/{kmers}/coverage.pmean.tsv"
    conda:
        "../../envs/bamm.yaml"
    shell: "bamm parse -c {output} -m pmean -b {input}"

rule mmgenome_orfs:
    input:
        "{project}/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.gz"
    output:
        orfs="{project}/mmgenome/{assembler}/{treatment}/{kmers}/orfs.faa",
        nucleotide="{project}/mmgenome/{assembler}/{treatment}/{kmers}/orfs.fna",
        orfscleaned="{project}/mmgenome/{assembler}/{treatment}/{kmers}/orfs.clean.faa"
    log:
        "{project}/mmgenome/{assembler}/{treatment}/{kmers}/prodigal.log"
    conda: 
        "../../envs/prodigal.yaml"
    shell: """
        zcat {input} | prodigal -d {output.nucleotide} -a {output.orfs} -i /dev/stdin -m -o {log} -p meta -q
        cut -f1 -d ' ' {output.orfs} > {output.orfscleaned}
    """

rule mmgenome_essential:
    input:
        "{project}/mmgenome/{assembler}/orfs.clean.faa"
    output:
        prediction="{project}/mmgenome/{assembler}/assembly.hmm.orfs.txt",
        essential="{project}/mmgenome/{assembler}/essential.txt"
    log:
       "{project}/mmgenome/{assembler}/prodigal.log"
    run:
        shell("hmmsearch --tblout {output.prediction} --cut_tc --notextw ~/install/mmgenome/scripts/essential.hmm {input} > {log}")
        shell("echo 'scaffold orf hmm.id' > {output.essential}")
        shell("tail -n+4 {output.prediction} | sed 's/ * / /g' | cut -f1,4 -d ' ' | sed 's/_/ /' >> essential.txt")

rule mmgenome_extract_essential:
    input:
        prediction="{project}/mmgenome/{assembler}/assembly.hmm.orfs.txt",
        orfs="{project}/mmgenome/{assembler}/orfs.clean.faa"
    output:
        posorfs="{project}/mmgenome/{assembler}/list.of.positive.orfs.txt",
        faa="{project}/mmgenome/{assembler}/assembly.orfs.hmm.faa"
    run:
        shell("grep -v '^#' {input.prediction} | cut -f1 -d ' ' > {output.posorfs}")
        shell("perl ~/install/mmgenome/scripts/extract.using.header.list.pl -l {output.posorfs} -s {input.orfs} -o {output.faa}")

# TODO: replace by diamond
# TODO: add the MEGAN command
rule mmgenome_essential_annotate:
    input:
        faa="{project}/mmgenome/{assembler}/assembly.orfs.hmm.faa"
    output:
        blast="{project}/mmgenome/{assembler}/assembly.orfs.hmm.blast.xml",
        tax="{project}/mmgenome/{assembler}/tax.txt"
    threads: 16
    run:
        shell("blastp -query {input.faa} -db /data/db/blast/nr/20150311/nr -evalue 1e-5 -num_threads {threads} -max_target_seqs 5 -outfmt 5 -out {output.blast}")
        # Here we need to run MEGAN first
        shell("perl ~/install/mmgenome/scripts/hmm.majority.vote.pl -i {output.blast} -o {output.tax}")

rule mmgenome_load_data:
     input:
         assembly="{project}/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.gz",
         essential="{project}/mmgenome/{assembler}/essential.txt",
         coverage="{project}/mmgenome/{assembler}/coverage.pmean.tsv",
         tax="{project}/mmgenome/{assembler}/tax.txt"
     output:
         "{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/{project}.RData"
     run:
         # TODO: make sample name independant
         R("""
library(mmgenome)
ess <- read.table("{input.essential}", header = T, sep = " ")
coverage = read.csv("{input.coverage}", header=T,sep="\t")
MPR1.1n.EST <- coverage[,c("X.contig","X1511KMI.0007.bamm.assembly.MPR1.1n.EST_1.bam")]
MPR1.1n.CUR <- coverage[,c("X.contig","X1511KMI.0007.bamm.assembly.MPR1.1n.CUR_1.bam")]
MPR2.1n.EST <- coverage[,c("X.contig","X1511KMI.0007.bamm.assembly.MPR2.1n.EST_1.bam")]

assembly <- readDNAStringSet("{input.assembly}", format = "fasta")
tax <- read.table("{input.tax}", header = T, sep = "\t")
d <- mmload(assembly = assembly, coverage = c("MPR1.1n.CUR", "MPR1.1n.EST", "MPR2.1n.EST"), essential=ess, tax=tax, tax.expand = "Proteobacteria", tax.freq = 85)
rm(list = c("assembly"))
save.image(file={output})
""")

