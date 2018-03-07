rule mmgenome_coverage:
    input:
        expand("{{project}}/bamm/{{assembler}}/{{treatment}}/{{kmers}}/assembly.{sample}_forward_paired.bam", sample=config["data"]) 
    output:
        "{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/coverage.pmean.tsv"
    conda:
        "../../envs/bamm.yaml"
    shell: "bamm parse -c {output} -m pmean -b {input}"

rule mmgenome_orfs:
    input:
        "{project}/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.gz"
    output:
        orfs="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/orfs.faa",
        nucleotide="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/orfs.fna",
        orfscleaned="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/orfs.clean.faa"
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
        "{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/orfs.clean.faa"
    output:
        prediction="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/assembly.hmm.orfs.txt",
        essential="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/essential.txt"
    log:
       "{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/prodigal.log"
    conda:
       "../../envs/mmgenome.yaml"
# TODO: remove hardcoded path
    shell: """
        hmmsearch --tblout {output.prediction} --cut_tc --notextw ~/install/mmgenome/scripts/essential.hmm {input} > {log}
        echo "scaffold orf hmm.id" > {output.essential}
        tail -n+4 {output.prediction} | sed "s/ * / /g" | cut -f1,4 -d " " | sed "s/_/ /" >> {output.essential}
        """

rule mmgenome_extract_essential:
    input:
        prediction="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/assembly.hmm.orfs.txt",
        orfs="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/orfs.clean.faa"
    output:
        posorfs="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/list.of.positive.orfs.txt",
        faa="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/assembly.orfs.hmm.faa"
    run:
        shell("grep -v '^#' {input.prediction} | cut -f1 -d ' ' > {output.posorfs}")
        shell("perl ~/install/mmgenome/scripts/extract.using.header.list.pl -l {output.posorfs} -s {input.orfs} -o {output.faa}")

# TODO: replace by diamond
# TODO: add the MEGAN command
# TODO: Update db ref and move to config
# https://madsalbertsen.github.io/multi-metagenome/docs/step5.html
rule mmgenome_essential_annotate:
    input:
        faa="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/assembly.orfs.hmm.faa"
    output:
        blast="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/assembly.orfs.hmm.blast.xml",
    threads: 16
    conda:
       "../../envs/mmgenome.yaml"
    message: "Finding essential genes - Extracting consensus taxonomic assignment"
    shell: """
        blastp -query {input.faa} -db /data/db/blast/nr/20150311/nr -evalue 1e-5 -num_threads {threads} -max_target_seqs 5 -outfmt 5 -out {output.blast}
        # Here we need to run MEGAN first
        java -Xmx32G -Djava.awt.headless=true -Duser.language=en -Duser.region=US -cp '/data/tools/MEGAN/6.10.8/jars/MEGAN.jar:/data/tools/MEGAN/6.10.8/jars/data.jar' megan.tools.Blast2LCA -i Shotgun_EPS/binning/mmgenome/megahit/all/longreads/assembly.orfs.hmm.blast.xml -f BlastXML -ms 50 -me 0.01 -top 50 -a2t /data/db/megan/prot_acc2tax-Oct2017X1.abin
        """
rule mmgenome_filter_megan:
    input:
        taxonomy="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/assembly.orfs.hmm.blast-taxonomy.txt",
    output:
        taxonomy="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/assembly.orfs.hmm.blast-taxonomy-filtered.txt"
    params:
        minscore="80"
    run:
        import re
        import sys
        
        minscore = int(params.minscore)
        out = open(output.taxonomy, 'w')

        filteredtax = {}
        # Parse the taxonomy string and apply score filter
        for line in open(input.taxonomy):
            taxdict = {}
            s = line.split(';')
            for i in range(2,len(s),2):                
                try:
                    level, value = s[i].strip().split('__')
                    if int(s[i+1]) >= 80:
                        taxdict[level] = s[i].strip()
                except:
                    pass
            taxonomy = ["root","cellular organisms"]
            for level in ["d","p","c","o","f","g","s"]:
                taxonomy.append(taxdict.setdefault(level, "unclassified"))

            taxstring = ";".join(taxonomy)
            genestring = ";".join(s[0:1])

            out.write("%s\t%s\n" % (genestring,taxstring))
        out.close()

rule mmgenome_consensus_tax:
   input:
       taxonomy="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/assembly.orfs.hmm.blast-taxonomy-filtered.txt"
   output:
       essential="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/tax.txt"
   shell: "perl ~/src/orochi/scripts/hmm.majority.vote.pl -i {input.taxonomy} -o {output.essential}"

rule mmgenome_load_data:
     input:
         assembly="{project}/assembly/{assembler}/{treatment}/{kmers}/assembly.fa.gz",
         essential="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/essential.txt",
         coverage="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/coverage.pmean.tsv",
         tax="{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/tax.txt"
     output:
         "{project}/binning/mmgenome/{assembler}/{treatment}/{kmers}/{project}.RData"
     conda:
         "../../envs/mmgenome.yaml"
     script:
         "mmgenome.R"
