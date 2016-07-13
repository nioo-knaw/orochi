__author__ = "Mattias de Hollander"
__copyright__ = "Copyright 2016, Mattias de Hollander"
__email__ = "m.dehollander@nioo.knaw.nl"
__license__ = "MIT"

from snakemake.utils import min_version
from snakemake.utils import R

min_version("3.5.4")

configfile: "config.json"

rule final:
        input: expand("{project}/trimming/{sample}_unpaired/forward_reads.fastq \
                       {project}/gunzip/{sample}.fastq \
                       {project}/stats/raw.readstat.csv \
                       {project}/stats/trimmed.readstat.csv \
                       {project}/diamond/{sample}.diamond.nr-taxonomy.tsv \
                       {project}/diamond/{sample}.diamond.nr-taxonomy-filtered.tsv \
                       {project}/diamond/{sample}.diamond.nr-taxonomy-filtered.qiime.txt \
                       {project}/diamond/{project}.RData \
                       {project}/extract_16S/{sample}.bbduk.fa.gz \
                       {project}/extract_16S/taxonomy/{project}.biom \
                       {project}/extract_16S/megahit/{sample}.gff \
                       {project}/extract_16S/megahit/{sample}/final.contigs.fa.rdp \
                       {project}/assembly/{assembler}/contigs.fasta \
                       {project}/stats/{assembler}.quast.report.txt \
                       {project}/stats/{assembler}.assembly.flagstat.txt \
                       {project}/microbecensus/{sample}.ags.txt \
                       {project}/genecatalog/all.coverage.norm.tpm.biom \
                       {project}/genecatalog/allgenecalled.diamond.nr-taxonomy.tsv \
                       {project}/genecatalog/all.coverage.norm.tpm.taxonomy.biom \
                       {project}/genecatalog/uproc/allgenecalled.uproc.kegg.txt \
                       {project}/genecatalog/uproc/allgenecalled.uproc.cog.txt \
                       {project}/genecatalog/uproc/allgenecalled.uproc.pfam.txt \
                       {project}/genecatalog/all.coverage.norm.tpm.taxonomy.ko.pfam.tsv".split(), project=config["project"], sample=config["data"], assembler=config["assembler"])

rule sickle_pe:
    input:
        forward = lambda wildcards: config["data"][wildcards.sample]['forward'],
        reverse = lambda wildcards: config["data"][wildcards.sample]['reverse']
    output:
        forward="{project}/trimming/{sample}_1.fastq.gz",
        reverse="{project}/trimming/{sample}_2.fastq.gz",
        unpaired="{project}/trimming/{sample}_unpaired.fastq.gz"
    params:
        dir = "sickle",
        quality=config['min_qual'],
        length=config['min_length']
    log:
       "{project}/trimming/{sample}.log"
    wrapper:
        "file://./bio/sickle_pe"

rule readstat_raw:
    input:
        expand("{{project}}/gunzip/{sample}.fastq", sample=config["data"])
    output:
        "{project}/stats/raw.readstat.csv"
    log:
        "{project}/stats/raw.readstat.log"
    shell: "set +u; source ~/.virtualenvs/khmer/bin/activate; set -u; /data/tools/khmer/scripts/readstats.py {input} --csv -o {output} 2> {log}"
 
rule readstat_trim:
    input:
        forward=expand("{{project}}/trimming/{sample}_1.fastq.gz", sample=config["data"]),
        reverse=expand("{{project}}/trimming/{sample}_2.fastq.gz", sample=config["data"]),
        unpaired=expand("{{project}}/trimming/{sample}_unpaired.fastq.gz", sample=config["data"]),
    output:
        "{project}/stats/trimmed.readstat.csv"
    log:
        "{project}/stats/trimmed.readstat.log"
    shell: "set +u; source ~/.virtualenvs/khmer/bin/activate; set -u; /data/tools/khmer/scripts/readstats.py {input} --csv -o {output} 2> {log}"

rule split_unpaired:
    input:
        "{project}/trimming/{sample}_unpaired.fastq.gz"
    output:
        "{project}/trimming/{sample}_unpaired/forward_reads.fastq",
        "{project}/trimming/{sample}_unpaired/reverse_reads.fastq"
    params:
        outdir = "{project}/trimming/{sample}_unpaired/"
    shell: "extract_reads_from_interleaved_file.py -i {input} -o {params.outdir}"

rule count_unpaired_forward:
    input:
        expand("{{project}}/trimming/{sample}_unpaired/forward_reads.fastq", sample=config["data"])
    output:
       "{project}/trimming/unpaired_forward.count.txt"
    params:
        samples=config["data"]
    run:
        for i, file in enumerate(input):
            sample = params.samples[i]
            shell("printf {sample}'\t' >> {output} && cat {file} | printf $((`wc -l`/4)) >> {output} && printf '\tforward\n' >> {output}")
            
rule count_unpaired_reverse:
    input:
        expand("{{project}}/trimming/{sample}_unpaired/reverse_reads.fastq", sample=config["data"])
    output:
       "{project}/trimming/unpaired_reverse.count.txt"
    params:
        samples=config["data"]
    run:
        for i, file in enumerate(input):
            sample = params.samples[i]
            shell("printf {sample}'\t' >> {output} && cat {file} | printf $((`wc -l`/4)) >> {output} && printf '\treverse\n' >> {output}")

rule diamond_per_sample:
    input:
        "{project}/trimming/{sample}_1.fastq.gz"
    output:
        tsv="{project}/diamond/{sample}.diamond.nr.tsv"
    params:
        reference=config["diamond_database"],
        version="0.7.10",
        output="{project}/diamond/{sample}.diamond.nr",
        format="tab",
        tmp="/tmp"        
    threads: 16
    run:
        shell("/data/tools/diamond/{params.version}/bin/diamond blastx -c 1 -d {params.reference} -t {params.tmp} -p {threads} -q {input} -a {params.output}")
        shell("/data/tools/diamond/{params.version}/bin/diamond view -f {params.format} -a {params.output}.daa -o {output.tsv}")

rule diamond_lca:
    input:
        "{project}/diamond/{sample}.diamond.nr.tsv"
    output:
        "{project}/diamond/{sample}.diamond.nr-taxonomy.tsv"
    # shell: "/data/tools/MEGAN/mtools/mtools/bin/lcamapper.sh -i {input} -f Detect -ms 50 -me 0.01 -tp 50 -gt /data/tools/MEGAN/mtools/mtools/data/gi_taxid_prot-4March2015.bin -o {output}"
    shell: "java -Xmx32G -Djava.awt.headless=true -Duser.language=en -Duser.region=US -cp '/data/tools/MEGAN/5.10.0/jars/MEGAN.jar:/data/tools/MEGAN/5.10.0/jars/data.jar' megan.tools.TaxonPathClassifier -i {input} -f Detect -ms 50 -me 0.01 -tp 50 -g2t /data/db/megan/gi_taxid_prot-4March2015.bin -o {output}"

rule seq_names_forward:
    input:
        "{project}/trimming/{sample}_1.fastq.gz"
    output:
        "{project}/trimming/{sample}_1.names.txt"
    shell: "zcat {input} | awk '(NR%4 == 1){{print $0}}' | cut -d' ' -f 1 | cut -c 2- > {output}" 

rule diamond_filter:
    input:
        taxonomy="{project}/diamond/{sample}.diamond.nr-taxonomy.tsv",
        reads = "{project}/trimming/{sample}_1.names.txt"
    output:
        taxonomy="{project}/diamond/{sample}.diamond.nr-taxonomy-filtered.tsv"
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
                    level,value = re.match('^\[(.*)\]\s(.*)', s[i].strip()).groups()                    
                    if int(s[i+1]) >= minscore:
                        taxdict[level] = value
                except:
                    pass
            taxonomy = ["Root"]
            # "Root", "k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae"
            #for level in ["SuperKingdom","Phylum","Class","Order","Family","Genus","Species"]:
            for level,prefix in [("SuperKingdom", "k__"),("Phylum", "p__"),("Class", "c__"),("Order", "o__"),("Family", "f__"),("Genus", "g__"),("Species", "s__")]:
                taxonomy.append(prefix + taxdict.setdefault(level, "unclassified"))

            taxstring = ";".join(taxonomy)
            genestring = ";".join(s[0:1])
            filteredtax[genestring] = taxstring

        #for gene in filteredtax:
        #    out.write("%s\t%s\n" % (gene,filteredtax[gene]))
        #out.close()

        # Add missing genes as unclassified
        genes = [line.strip().split()[0] for line in open(input.reads).readlines()][1:]        
        for gene in genes:
            if gene not in filteredtax:
                filteredtax[gene] = 'Root;k__unclassified;p__unclassified;c__unclassified;o__unclassified;f__unclassified;g__unclassified;s__unclassified;'        
            out.write("%s\t%s\n" % (gene,filteredtax[gene]))
        out.close()

rule diamond2phyloseq:
    input:
        "{project}/diamond/{sample}.diamond.nr-taxonomy-filtered.tsv"
    output:
        "{project}/diamond/{sample}.diamond.nr-taxonomy-filtered.qiime.txt"
    params:
        sample = "{sample}"
    run:
        shell("echo -e \"#OTU ID\t{params.sample}\tConsensus Lineage\" > {output}")
        # Add a column with count is one, because lines are raw reads (e.g. no clustering)
        shell("awk 'BEGIN {{OFS=\"\t\"}} {{ print $1,1,$2 }}' {input} >> {output}")

rule diamond_R:
    input:
        expand("{{project}}/diamond/{sample}.diamond.nr-taxonomy-filtered.qiime.txt", sample=config["data"])
    output:
        "{project}/diamond/{project}.RData"
    run:
        files = "','".join(input)
        vector = "c('%s')" % files

        filelist = []
        for i in range(1,len(input)+1):
            filelist.append("d[[%i]]"  % i)
        filestring = ",".join(filelist)
        print(filestring)

        R("""
        library(phyloseq)
        library(parallel)

        f <- eval(parse(text="{vector}"))
        d <- mclapply(f, import_qiime, mc.cores=18)
        #p <- merge_phyloseq(d[[6]],d[[7]],d[[8]],d[[9]],d[[10]],d[[11]],d[[12]],d[[13]],d[[14]],d[[15]],d[[16]],d[[17]],d[[18]])
        p <- merge_phyloseq({filestring})

        glom_kingdom <- tax_glom(p, taxrank = 'Kingdom')
        glom_phylum <- tax_glom(p, taxrank = 'Phylum')
        glom_class <- tax_glom(p, taxrank = 'Class') 
        glom_order <- tax_glom(p, taxrank = 'Order')
        glom_family <- tax_glom(p, taxrank = 'Family')
        glom_genus <- tax_glom(p, taxrank = 'Genus')
        glom_species <- tax_glom(p, taxrank = 'Species')


        rm(list = c("f", "d", "p"))
        save.image(file="{output}")
        """)

rule interleave:
    input:
#        "{project}/trimming/{sample}_1.fastq.gz",
#        "{project}/trimming/{sample}_2.fastq.gz"
        forward = lambda wildcards: config["data"][wildcards.sample]['forward'],
        reverse = lambda wildcards: config["data"][wildcards.sample]['reverse']
    output:
        "{project}/gunzip/{sample}.fastq"
    shell: "set +u; source ~/.virtualenvs/khmer/bin/activate; set -u; /data/tools/khmer/scripts/interleave-reads.py -o {output} {input}"

rule interleave_merge:
    input:
        fasta = expand("{{project}}/interleaved/{sample}.fastq", sample=config["data"]),
    output: 
        fasta="{project}/interleaved_merged/all.fastq"
    shell: """cat {input} > {output}"""


rule kraken_per_sample:
    input:
        "{project}/trimming/{sample}_1.fastq.gz",
        "{project}/trimming/{sample}_2.fastq.gz"
    output:
        kraken = "{project}/kraken/{sample}.kraken",
        taxonomy = "{project}/kraken/{sample}.kraken.taxonomy",
        report = "{project}/kraken/{sample}.kraken.report"
    log:
        "{project}/kraken/{sample}.kraken.log"
    threads: 16
    run:
        shell("/data/tools/kraken/0.10.5-beta/bin/kraken --db /data/db/kraken/ --threads {threads} --fastq-input {input} --paired --check-names --gzip-compressed --output {output.kraken} 2> {log}")
        shell("/data/tools/kraken/0.10.5-beta/bin/kraken-translate --db /data/db/kraken/ --mpa-format {output.kraken} > {output.taxonomy}")
        shell("/data/tools/kraken/0.10.5-beta/bin/kraken-mpa-report --db /data/db/kraken/  {output.kraken} > {output.report}")

#TODO:
# Make plots for raw reads tax assignment with kraken/diamond
# Convert to a dummy otutable to load with phyloseq with a fake count column
# Add qiime header (todo)
# awk 'BEGIN {OFS="\t"} { print $1,1,$2 }' 1511KMI-0007/kraken/MPR1-1n-EST.kraken.taxonomy | sed 's/|/;/g' > test2.qiime
# make plots with gpplot

rule kraken2phyloseq:
    input:
        "{project}/kraken/{sample}.kraken.taxonomy"
    output:
        "{project}/kraken/{sample}.kraken.qiime.taxonomy"
    params:
        sample = "{sample}"
    run:
        shell("echo -e \"#OTU ID\t{params.sample}\tConsensus Lineage\" > {output}")
        # Add a column with count is one, because lines are raw reads (e.g. no clustering)
        # Replace pipe sign with semicolon for qiime compatebility 
        shell("awk 'BEGIN {{OFS=\"\t\"}} {{ print $1,1,$2 }}' {input} | sed 's/|/;/g' >> {output}")

rule extract_16S:
    input:
        forward="{project}/trimming/{sample}_1.fastq.gz",
        reverse="{project}/trimming/{sample}_2.fastq.gz"
    output:
        "{project}/extract_16S/{sample}.bbduk.fa.gz"
    log:
        "{project}/extract_16S/{sample}.bbduk.log"
    threads: 16
    shell: "/data/tools/bbmap/34.08/bbduk.sh in={input.forward} in2={input.reverse} outm={output} k=31 ref=/data/db/Silva/ribokmers.fa.gz threads={threads} 2> {log}"

rule rename_16S:
    input:
        "{project}/extract_16S/{sample}.bbduk.fa.gz"
    output:
        "{project}/extract_16S/{sample}.bbduk.rename.fa.gz"
    params:
        prefix="{sample}"
    threads: 2
    run: 
        # Quick fix because uc2otuable does not allow underscore in sample names
        sample = params.prefix.replace("_","-")
        shell("zcat {input} | awk '/^>/ {{$0=\">{sample}_\" substr($0,2)}}1' | gzip > {output}")

rule merge_16S:
    input:
        fasta = expand("{{project}}/extract_16S/{sample}.bbduk.rename.fa.gz", sample=config["data"]),
    output: 
        fasta="{project}/extract_16S/all.fa.gz"
    shell: """zcat {input} | gzip > {output}"""

rule map_to_16S:
    input:
        "{project}/extract_16S/all.fa.gz"
    output:
        mapping="{project}/extract_16S/mapping/{project}.uc",
        biom="{project}/extract_16S/otutable/{project}.biom",
        biom_tax="{project}/extract_16S/taxonomy/{project}.biom"
    log:
        "{project}/extract_16S/mapping/vsearch.log"
    params:
        silva_ref=config['silva_ref'],
        silva_taxonomy=config['silva_taxonomy']
    threads: 16
    run:
        shell("source /data/tools/vsearch/1.9.6/env.sh; vsearch --usearch_global {input} -db {params.silva_ref} -id 0.97 --strand both --uc {output.mapping} --threads {threads} --log {log}")
        shell("set +u; source /data/tools/biom-format/2.1.5/env.sh; set -u; biom from-uc -i {output.mapping} -o {output.biom}")
        shell("set +u; source /data/tools/qiime/1.9.1/env.sh; set -u; biom add-metadata -i {output.biom} -o {output.biom_tax} --observation-metadata-fp {params.silva_taxonomy} --observation-header OTUID,taxonomy --sc-separated taxonomy --float-fields confidenc")

rule megahit_16S:
    input:
        "{project}/extract_16S/{sample}.bbduk.fa.gz"
    output:
        contigs="{project}/extract_16S/megahit/{sample}/final.contigs.fa"
    params: 
        dir="{project}/extract_16S/megahit/{sample}/"
    threads: 8
    shell: "/data/tools/megahit/1.0.3/megahit --continue -t {threads} --12 {input} -o {params.dir}"

rule barrnap:
    input:
        "{project}/extract_16S/megahit/{sample}/final.contigs.fa"
    output:
        "{project}/extract_16S/megahit/{sample}.gff"
    log: "{project}/extract_16S/megahit/{sample}/barrnap.log"
    threads: 16
    run:
        shell("/data/tools/barrnap/0.6/bin/barrnap --threads {threads} --reject 0.1 {input} > {output} 2> {log}")

rule rdp_16S:
    input:
        contigs="{project}/extract_16S/megahit/{sample}/final.contigs.fa"
    output:
        contigs="{project}/extract_16S/megahit/{sample}/final.contigs.fa.rdp"
    shell: "java -Xmx1g -jar /data/tools/rdp-classifier/2.10/classifier.jar classify -c 0.8 {input} -f filterbyconf -o {output}"

rule megahit_16S_cross_assembly:
    input:
        expand("{{project}}/extract_16S/{sample}.bbduk.fa.gz", sample=config["data"])
    output:
        contigs="{project}/extract_16S/megahit/cross_assembly/final.contigs.fa"
    params:
        dir="{project}/extract_16S/megahit/cross_assembly/"
    threads: 8
    run:
        input_str = ",".join(input) 
        shell("/data/tools/megahit/1.0.3/megahit --continue -t {threads} --12 {input_str} -o {params.dir}")

rule barrnap_cross_assembly:
    input:
        "{project}/extract_16S/megahit/cross_assembly/final.contigs.fa"
    output:
        "{project}/16S_assembly/barrnap.gff"
    log: "{project}/16S_assembly/barrnap.log"
    threads: 16
    run:
        shell("/data/tools/barrnap/0.7_gene/barrnap/bin/barrnap --gene ssu --threads {threads} --reject 0.1 {input} > {output} 2> {log}")

rule rdp_16S_cross_assembly:
    input:
        contigs="{project}/extract_16S/megahit/cross_assebmly/final.contigs.fa"
    output:
        contigs="{project}/extract_16S/megahit/cross_assembly/final.contigs.fa.rdp"
    shell: "java -Xmx1g -jar /data/tools/rdp-classifier/2.10/classifier.jar classify -c 0.8 {input} -f filterbyconf -o {output}"

rule combine_reads:
    input:
        forward=expand("{{project}}/trimming/{sample}_1.fastq.gz", sample=config["data"]),
        reverse=expand("{{project}}/trimming/{sample}_2.fastq.gz", sample=config["data"]),
        unpaired=expand("{{project}}/trimming/{sample}_unpaired.fastq.gz", sample=config["data"])
    output:
        forward="{project}/trimming/R1.fastq.gz",
        reverse="{project}/trimming/R2.fastq.gz",
        unpaired="{project}/trimming/singular.fastq.gz"
    run:
        shell("zcat `echo {input.forward} | xargs -n1 | sort | xargs` | gzip > {output.forward}")
        shell("zcat `echo {input.reverse} | xargs -n1 | sort | xargs` | gzip > {output.reverse}")
        shell("zcat `echo {input.unpaired} | xargs -n1 | sort | xargs` | gzip > {output.unpaired}")

rule megahit:
    input:
        forward=expand("{{project}}/trimming/{sample}_1.fastq.gz", sample=config["data"]),
        reverse=expand("{{project}}/trimming/{sample}_2.fastq.gz", sample=config["data"]),
        unpaired=expand("{{project}}/trimming/{sample}_unpaired.fastq.gz", sample=config["data"])
    output:
        contigs="{project}/megahit/final.contigs.fa",
        contigs_gzip="{project}/megahit/final.contigs.fa.gz",
        # This file contains all the settings of a run. When this file is not present megahit with run in normal mode, otherwise it continues with previous settings
        opts="{project}/megahit/opts.txt"
    params: dir="{project}/megahit"
    log: "{project}/megahit/megahit.log"
    threads: 16
    run:
        forward_str = ",".join(input.forward)
        reverse_str = ",".join(input.reverse) 
        unpaired_str = ",".join(input.unpaired) 

        # Parameter settings
        # meta            '--min-count 2 --k-list 21,41,61,81,99'             (generic metagenomes, default)
        # meta-sensitive  '--min-count 2 --k-list 21,31,41,51,61,71,81,91,99' (more sensitive but slower)
        # meta-large      '--min-count 2 --k-list 27,37,47,57,67,77,87'       (large & complex metagenomes, like soil)
        # bulk            '--min-count 3 --k-list 31,51,71,91,99 --no-mercy'  (experimental, standard bulk sequencing with >= 30x depth)
        # single-cell     '--min-count 3 --k-list 21,33,55,77,99,121 --merge_level 20,0.96' (experimental, single cell data)

        shell("/data/tools/megahit/1.0.6/megahit --continue --out-dir {params.dir} -m 0.9 --max-read-len 302 --cpu-only -t {threads} --presets meta -1 {forward_str} -2 {reverse_str} -r {unpaired_str} 2> {log}")
        shell("gzip -c {output.contigs} > {output.contigs_gzip}")

rule spades:
   input:
        forward=expand("{{project}}/trimming/{sample}_1.fastq.gz", sample=config["data"]),
        reverse=expand("{{project}}/trimming/{sample}_2.fastq.gz", sample=config["data"]),
        unpaired=expand("{{project}}/trimming/{sample}_unpaired.fastq.gz", sample=config["data"])
   output:
        "{project}/assembly/spades/contigs.fasta"
   params:
       outdir="{project}/assembly/spades/"
   log:
       "{project}/assembly/spades/spades.log"
   threads: 16
   run:
       forward_str = " -1 ".join(input.forward)
       reverse_str = " -2 ".join(input.reverse) 
       unpaired_str = " -s ".join(input.unpaired)  
       shell("ulimit -m 480000000; /data/tools/SPAdes/3.8.2/bin/spades.py -m 480 -1 {forward_str} -2 {reverse_str} -s {unpaired_str} --only-assembler -t {threads} -k 21,33,55,77 --careful -o {params.outdir} --tmp-dir {params.outdir}/tmp/ 2>&1 > /dev/null")

# Interleave paired end reads and convert to fasta
rule idba_prepare:
    input:
        forward="{project}/trimming/{sample}_1.fastq.gz",
        reverse="{project}/trimming/{sample}_2.fastq.gz",   
    output:
        forward=temp("{project}/trimming/{sample}_1.fastq"),
        reverse=temp("{project}/trimming/{sample}_2.fastq"),
        merged="{project}/trimming/{sample}.fasta"
    run:
        shell("gunzip -c {input.forward} > {output.forward}")
        shell("gunzip -c {input.reverse} > {output.reverse}")
        shell("/data/tools/idba/1.1.3/bin/fq2fa --merge {output.forward} {output.reverse} {output.merged}")

rule idba:
    input:
        expand("{{project}}/trimming/{sample}.fasta", sample=config["data"])
    output:
        "{project}/assembly/idba/scaffold.fa"
    params:
        outdir="{project}/assembly/idba/"
    threads: 25
    shell: "cat {input} | /data/tools/idba/1.1.3/bin/idba_ud -r /dev/stdin -o {params.outdir} --num_threads {threads} --mink 21 --maxk 21 --step 20 --pre_correction"

rule quast:
    input:
        "{project}/assembly/{assembler}/assembly.fa.gz"
    output:
        quast="{project}/assembly/{assembler}/quast/report.txt",
        stats="{project}/stats/{assembler}.quast.report.txt"
    params:
        outdir="{project}/assembly/{assembler}/quast"
    log:
        "{project}/assembly/{assembler}/quast/quast.log"
    threads: 16
    run:
        shell("python2.7 /data/tools/quast/4.1/bin/metaquast.py -o {params.outdir} --min-contig 0 --max-ref-number 0 -t {threads} {input} 2>&1 > {log}")
        shell("cp {output.quast} {output.stats}")

rule barrnap_cross_assembly_all:
    input:
        "{project}/megahit/assembly.fa"
    output:
        gff="{project}/megahit/barrnap.gff",
        fasta="{project}/megahit/barrnap.fasta"
    log: "{project}/megahit/barrnap.log"
    threads: 16
    run:
        shell("/data/tools/barrnap/0.7_gene/barrnap/bin/barrnap --gene ssu --threads {threads} --reject 0.1 {input} > {output.gff} 2> {log}")
        shell("bedtools getfasta -fi {input} -bed {output.gff} -fo {output.fasta}")

rule bamm:
    input:
        contigs="{project}/megahit/final.contigs.fa.gz",
        forward="{project}/trimming/{sample}_1.fastq.gz",
        reverse="{project}/trimming/{sample}_2.fastq.gz",
        unpaired="{project}/trimming/{sample}_unpaired.fastq.gz",
        index="{project}/megahit/final.contigs.fa.gz.bwt"
    output:
        "{project}/bamm/final.contigs.{sample}_1.bam",
        "{project}/bamm/final.contigs.{sample}_1.bam.bai"
    log:
        "{project}/bamm/{sample}.log"
    params:
        outdir="{project}/bamm/"
    threads: 16
    shell: "source /data/tools/samtools/1.3/env.sh; source /data/tools/BamM/1.7.0/env.sh; bamm make --kept -d {input.contigs} -c {input.forward} {input.reverse} -s {input.unpaired} -o {params.outdir} -t {threads} 2> {log}"

# This rule also produces bam files per sample... Not desired result.
rule bamm_all:
    input:
        contigs="{project}/megahit/final.contigs.fa.gz",
        forward=expand("{{project}}/trimming/{sample}_1.fastq.gz", sample=config["data"]),
        reverse=expand("{{project}}/trimming/{sample}_2.fastq.gz",  sample=config["data"]),
        unpaired=expand("{{project}}/trimming/{sample}_unpaired.fastq.gz",  sample=config["data"])
    output:
        "{project}/bammall/final.contigs.all.bam",
        "{project}/bammall/final.contigs.all.bam.bai"
    log:
        "{project}/bamm/all.log"
    params:
        outdir="{project}/bamm/"
    threads: 16
    run: 
        #make pairs of forward and reverse fastq files
        for i, element in enumerate(input.reverse, 1):
            input.forward.insert(i * 2 - 1, element)
        shell("source /data/tools/BamM/dev/env.sh; bamm make --force -d {input.contigs} -c {input.forward} -s {input.unpaired} -o {params.outdir} -t {threads} 2> {log}")

# samtools merge all.bam *.bam
# bamm parse -c coverage.tsv -m counts -b 1511KMI-0007/bamm/all.bam -t 16

rule bwa_index:
    input:
         "{project}/megahit/final.contigs.fa.gz"
    output:
        "{project}/megahit/final.contigs.fa.gz.bwt"
    log: "{project}/megahit/bwa-index.log"
    shell: "/data/tools/bwa/default/bin/bwa index {input} > {log}"


#
# Use the bwa and samtools wrapper from the snakemake-wrappers repository
#
# TODO: Unpaired reads can not be mapped at the same time
rule bwa_mem:
    input:
        ref="{project}/megahit/final.contigs.fa.gz",
        index="{project}/megahit/final.contigs.fa.gz.bwt",
        forward="{project}/trimming/{sample}_1.fastq.gz",
        reverse="{project}/trimming/{sample}_2.fastq.gz",
        unpaired="{project}/trimming/{sample}_unpaired.fastq.gz",
        sample=["{project}/trimming/{sample}_1.fastq.gz", "{project}/trimming/{sample}_2.fastq.gz"]
    output:
        "{project}/mapped/{sample}.bam"
    log:
        "{project}/logs/bwa_mem/{sample}.log"
#    params:
#        "-R '@RG\tID:{sample}\tSM:{sample}'"  # optional parameters for bwa mem (e.g. read group)
    threads: 8
    wrapper:
        "0.0.11/bio/bwa_mem"

rule bwa_mem_unpaired:
    input:
        ref="{project}/megahit/final.contigs.fa.gz",
        index="{project}/megahit/final.contigs.fa.gz.bwt",
        sample="{project}/trimming/{sample}_unpaired.fastq.gz",
    output:
        "{project}/mapped/{sample}_unpaired.bam"
    log:
        "{project}/logs/bwa_mem/{sample}.log"
    threads: 8
    wrapper:
        "0.0.11/bio/bwa_mem"

rule samtools merge:
    input:
        expand("{{project}}/bamm/{{assembler}}/assembly.{sample}_1.bam", sample=config["data"])
    output:
        "{project}/bamm/{assembler}/assembly.bam"
    shell:
        "samtools merge {output} {input}"
    
rule samtools_flagstat:
    input:
        "{project}/bamm/{assembler}/assembly.bam"
    output:
        "{project}/stats/{assembler}.assembly.flagstat.txt"
    shell:
        "samtools flagstat {input} > {output}"

#
# mmgenome
#

rule prepare_mmgenome:
    input:
        "{project}/assembly/megahit/final.contigs.fa.gz"
    output:
        gzip="{project}/assembly/megahit/assembly.fa.gz",
        fasta=temp("{project}/assembly/megahit/assembly.fa")
    run:
        shell("zcat {input} | awk '{{print $1}}' | sed 's/_/contig/' > {output.fasta}")
        shell("gzip -c {output.fasta} > {output.gzip}")

rule prepare_mmgenome_spades:
    input:
        "{project}/assembly/spades/contigs.fasta"
    output:
        gzip="{project}/assembly/spades/assembly.fa.gz",
        fasta=temp("{project}/assembly/spades/assembly.fa")
    run:
        shell("cat {input} | awk '{{print $1}}' | sed 's/_/contig/' > {output.fasta}")
        shell("gzip -c {output.fasta} > {output.gzip}")


rule mmgenome_bwa_index:
    input:
         "{project}/assembly/{assembler}/assembly.fa.gz"
    output:
        "{project}/assembly/{assembler}/assembly.fa.gz.bwt"
    log: "{project}/assembly/{assembler}/bwa-index.log"
    shell: "/data/tools/bwa/default/bin/bwa index {input} > {log}"

rule bamm_mmgenome:
    input:
        contigs="{project}/assembly/{assembler}/assembly.fa.gz",
        forward="{project}/trimming/{sample}_1.fastq.gz",
        reverse="{project}/trimming/{sample}_2.fastq.gz",
        unpaired="{project}/trimming/{sample}_unpaired.fastq.gz",
        index="{project}/assembly/{assembler}/assembly.fa.gz.bwt"
    output:
        "{project}/bamm/{assembler}/assembly.{sample}_1.bam",
        "{project}/bamm/{assembler}/assembly.{sample}_1.bam.bai"
    log:
        "{project}/bamm/{assembler}/{sample}.log"
    params:
        outdir="{project}/bamm/{assembler}"
    threads: 16
    shell: "set +u ;source /data/tools/samtools/1.3/env.sh; source /data/tools/BamM/1.7.0/env.sh; set -u; bamm make --keep_unmapped --kept -d {input.contigs} -c {input.forward} {input.reverse} -s {input.unpaired} -o {params.outdir} -t {threads} 2> {log}"

rule mmgenome_coverage:
    input:
        expand("{{project}}/bamm/assembly.{sample}_1.bam", sample=config["data"])
    output:
        "{project}/mmgenome/coverage.pmean.tsv"
    shell: "source /data/tools/samtools/1.3/env.sh; source /data/tools/BamM/1.7.0/env.sh; bamm parse -c {output} -m pmean -b {input}"

rule mmgenome_orfs:
    input:
        "{project}/assembly/{assembler}/assembly.fa.gz"
    output:
        orfs="{project}/mmgenome/{assembler}/orfs.faa",
        nucleotide="{project}/mmgenome/{assembler}/orfs.fna",
        orfscleaned="{project}/mmgenome/{assembler}/orfs.clean.faa"
    log:
        "{project}/mmgenome/{assembler}/prodigal.log"
    run:
        shell("zcat {input} | prodigal -d {output.nucleotide} -a {output.orfs} -i /dev/stdin -m -o {log} -p meta -q")
        shell("cut -f1 -d ' ' {output.orfs} > {output.orfscleaned}")

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
         assembly="{project}/assembly/{assembler}/assembly.fa.gz",
         essential="{project}/mmgenome/{assembler}/essential.txt",
         coverage="{project}/mmgenome/{assembler}/coverage.pmean.tsv",
         tax="{project}/mmgenome/{assembler}/tax.txt"
     output:
         "{project}/mmgenome/{project}.RData"
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

rule metabat:
    input: 
        contigs="{project}/assembly/{assembler}/assembly.fa",
        #bam=expand("{{project}}/mapped/{sample}.bam", sample=config["data"]),
        bam=expand("{{project}}/bamm/{{assembler}}/assembly.{sample}_1.bam", sample=config["data"])
    output:
        depth="{project}/metabat/depth.txt",
        bin="{project}/metabat/bin.1.fa"
    params:
        prefix="{project}/metabat/bin"
    log: "{project}/metabat/metabat.log"
    threads: 16
    run: 
        shell("/data/tools/metabat/dev/bin/jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}")
        shell("/data/tools/metabat/dev/bin/metabat -i {input.contigs} -a {output.depth} -o {params.prefix} --very-sensitive --numThreads {threads} --minContigByCorr 1500 --saveTNF saved.tnf --saveDistance saved.dist --verbose > {log}")

rule mmgenome_metabat:
    input:
        bin="{project}/metabat/bin.1.fa"
    output:
        "{project}/mmgenome/metabat.bins.txt"
    params:
        dir="{project}/metabat/"
    run:
        shell("echo -e 'scaffold\tbin' > {output}")
        shell("for file in `ls {params.dir}/*.fa` ; do noext=${{file%.fa}}; bin=$(basename ${{noext}}); awk -v bin=$bin '/^>/ {{split(substr($0,2),a,\":\"); print a[1] \"\\t\" bin;}}' $file;  done >> {output}")
        shell("sed --in-place -e s/[.]//g {output}")

rule checkm_lineage_metabat:
    input:
        bin="{project}/metabat/bin.1.fa"
    output:
        log="{project}/metabat/checkm/log.txt"
    params:
        indir="{project}/metabat/",
        outdir="{project}/metabat/checkm/"
    threads: 16
    shell: "set +u; source ~/.virtualenvs/groopm/bin/activate; source /data/tools/CheckM/0.9.7/env.sh; source /data/tools/pplacer/1.1/env.sh; set -u; python /data/tools/CheckM/0.9.7/bin/checkm lineage_wf -x fa -t {threads} -f {output.log} {params.indir} {params.outdir}"

rule megahit_per_sample:
    input:
        forward="{project}/trimming/{sample}_1.fastq.gz",
        reverse="{project}/trimming/{sample}_2.fastq.gz",
        unpaired="{project}/trimming/{sample}_unpaired.fastq.gz"
    output:
        contigs="{project}/megahit_per_sample/{sample}/final.contigs.fa",
        contigs_gzip="{project}/megahit_per_sample/{sample}/final.contigs.fa.gz",
        # This file contains all the settings of a run. When this file is not present megahit with run in normal mode, otherwise it continues with previous settings
        opts="{project}/megahit_per_sample/{sample}/opts.txt"
    params: dir="{project}/megahit_per_sample/{sample}/"
    log: "{project}/megahit_per_sample/{sample}/megahit.log"
    threads: 16
    run:
        # Parameter settings
        # meta            '--min-count 2 --k-list 21,41,61,81,99'             (generic metagenomes, default)
        # meta-sensitive  '--min-count 2 --k-list 21,31,41,51,61,71,81,91,99' (more sensitive but slower)
        # meta-large      '--min-count 2 --k-list 27,37,47,57,67,77,87'       (large & complex metagenomes, like soil)
        # bulk            '--min-count 3 --k-list 31,51,71,91,99 --no-mercy'  (experimental, standard bulk sequencing with >= 30x depth)
        # single-cell     '--min-count 3 --k-list 21,33,55,77,99,121 --merge_level 20,0.96' (experimental, single cell data)

        shell("/data/tools/megahit/1.0.3/megahit --continue --out-dir {params.dir} -m 0.9 --max-read-len 302 --cpu-only -t {threads} --presets meta -1 {input.forward} -2 {input.reverse} -r {input.unpaired} 2> {log}")
        shell("gzip -c {output.contigs} > {output.contigs_gzip}")

rule barrnap_per_sample:
    input:
        "{project}/megahit_per_sample/{sample}/final.contigs.fa"
    output:
        "{project}/megahit_per_sample/{sample}/barrnap.gff"
    log: "{project}/megahit_per_sample/{sample}/barrnap.log"
    threads: 16
    run:
        shell("/data/tools/barrnap/0.6/bin/barrnap --threads {threads} --reject 0.1 {input} > {output} 2> {log}")

rule bamm_per_sample:
    input:
        contigs="{project}/megahit_per_sample/{sample}/final.contigs.fa.gz",
        forward="{project}/trimming/{sample}_1.fastq.gz",
        reverse="{project}/trimming/{sample}_2.fastq.gz",
        unpaired="{project}/trimming/{sample}_unpaired.fastq.gz"
    output:
        "{project}/bamm/{sample}/final.contigs.{sample}_1.bam",
        "{project}/bamm/{sample}/final.contigs.{sample}_1.bam.bai"
    log:
        "{project}/bamm/{sample}.log"
    params:
        outdir="{project}/bamm/{sample}"
    threads: 16
    shell: "source /data/tools/samtools/1.3/env.sh; source /data/tools/BamM/1.7.0/env.sh; bamm make --force -d {input.contigs} -c {input.forward} {input.reverse} -s {input.unpaired} -o {params.outdir} --keep_unmapped -t {threads} 2> {log}"

# zcat 1511KMI-0007/megahit_per_sample/MPR1-1n-CUR/final.contigs.fa.gz | gmhmmp -a -d -f G -m /data/tools/metagenemark/3.26/model/MetaGeneMark_v1.mod -o test.gff /dev/stdin
# zcat 1511KMI-0007/megahit_per_sample/MPR1-1n-CUR/final.contigs.fa.gz | /data/tools/prodigal/2.62/bin/prodigal -d prodigal.fna -a prodigal.faa -i /dev/stdin -m -f gff -o prodigal.gff -p meta

rule genemark:
    input:
        #contigs="{project}/megahit_per_sample/{sample}/final.contigs.fa.gz",
        contigs="{project}/megahit/final.contigs.fa.gz"
    output:
        gff="{project}/genemark/{sample}.gff",
        nucleotide="{project}/genemark/{sample}.fna",
        protein="{project}/genemark/{sample}.faa"
    shell: "source /data/tools/metagenemark/3.26/env.sh; zcat {input} | gmhmmp -A {output.protein} -D {output.nucleotide} -f G -m /data/tools/metagenemark/3.26/model/MetaGeneMark_v1.mod -o {output.gff} /dev/stdin"

#TODO: Fix header, few duplicates present after merging    
rule orfs_unique:
    input:
        nucleotide="{project}/mmgenome/orfs.fna",
        protein="{project}/mmgenome/orfs.faa",
    output:
        nucleotide="{project}/mmgenome/orfs.fna.gz",
        protein="{project}/mmgenome/orfs.faa.gz"
    params:
        # TODO: does this needs to be changed?
        prefix="prefix"
    run:
        # Headers of the nuc and aa files need to be the same in order to construct the an aa file from the nucl genecatalog  
        # add a prefix to the header, keep only the part until the first |, remove empty lines with sed
        shell("cat {input.nucleotide} | awk -F '|' '/^>/ {{$0=\">{params.prefix}_\" substr($1,2)}}1' | sed '/^\s*$/d' |  gzip > {output.nucleotide}")
        shell("cat {input.protein} | awk -F '|' '/^>/ {{$0=\">{params.prefix}_\" substr($1,2)}}1' | sed '/^\s*$/d' | gzip > {output.protein}")


rule genemark_merge:
    input:
        nucleotide=expand("{{project}}/genemark/{sample}.fna.gz", sample=config["data"]),
        protein=expand("{{project}}/genemark/{sample}.faa.gz", sample=config["data"])
    output:
        # Place merged result in a new directory otherwise it matches with the sample names
        nucleotide="{project}/genemark_merged/all.fna.gz",
        protein="{project}/genemark_merged/all.faa.gz"
    run:
        shell("zcat {input.nucleotide} | gzip > {output.nucleotide}")
        shell("zcat {input.protein} | gzip > {output.protein}")

rule genecatalog:
    input:
        "{project}/mmgenome/orfs.fna.gz"
    output:
        clusters="{project}/genecatalog/allgenecalled.clusters.uc",
        fasta="{project}/genecatalog/allgenecalled.centroids.fna"
    threads: 16
    shell: "vsearch -cluster_fast {input} -id 0.95 -idprefix 4 -target_cov .9 -consout {output.fasta} -uc {output.clusters} -threads {threads}"

rule genecatalog_aa:
    input:
        genecatalog="{project}/genecatalog/allgenecalled.centroids.fna",
        protein="{project}/mmgenome/orfs.faa.gz"
    output:
        nucleotide="{project}/genecatalog/allgenecalled.fna.gz",
        protein="{project}/genecatalog/allgenecalled.faa.gz"
    run: 
        # Select the headers of the nucleotide sequences in the genecatalog
        # Select the seqids inbetween the = and ;
        # Select all the protein sequences with the same seqids
        shell("grep '>' {input.genecatalog} | awk -F'[=;]' '{{print $2}}' | seqtk subseq {input.protein} /dev/stdin | gzip > {output.protein}")
        # Make sure the nucleotide sequences have the same as the proteins
        shell("awk -F'[=;]' '/^>/ {{$0=\">\"$2}}1' {input.genecatalog} | gzip > {output.nucleotide}")

# bwa index 1511KMI-0007/genecatalog/allgenecalled.centroids.fa

rule genecatalog_bwa_index:
   input:
       genes="{project}/genecatalog/allgenecalled.fna.gz"
   output:
       index="{project}/genecatalog/allgenecalled.fna.gz.bwt"
   shell: "bwa index {input}"

rule map_to_genes:
    input:
        genes="{project}/genecatalog/allgenecalled.fna.gz",
        index="{project}/genecatalog/allgenecalled.fna.gz.bwt",
        forward="{project}/trimming/{sample}_1.fastq.gz",
        reverse="{project}/trimming/{sample}_2.fastq.gz",
        unpaired="{project}/trimming/{sample}_unpaired.fastq.gz"
    output:
        "{project}/genecatalog/{sample}/allgenecalled.{sample}_1.bam",
        "{project}/genecatalog/{sample}/allgenecalled.{sample}_1.bam.bai",
        "{project}/genecatalog/{sample}/allgenecalled.{sample}_unpaired.bam",
        "{project}/genecatalog/{sample}/allgenecalled.{sample}_unpaired.bam.bai"

    log:
        "{project}/genecatalog/{sample}/mapping.log"
    params:
        outdir="{project}/genecatalog/{sample}"
    threads: 16
    # Use --kept to re/multi use preindexed reference. Otherwise with --force the indexes are rebuild every time
    shell: "source /data/tools/samtools/1.3/env.sh; source /data/tools/BamM/1.7.0/env.sh; bamm make --kept -d {input.genes} -c {input.forward} {input.reverse} -s {input.unpaired} -o {params.outdir} --keep_unmapped -t {threads} 2> {log}"

rule coveragetable: 
    input: 
        paired = "{project}/genecatalog/{sample}/allgenecalled.{sample}_1.bam",
        unpaired = "{project}/genecatalog/{sample}/allgenecalled.{sample}_unpaired.bam"
    output:
        paired = "{project}/genecatalog/{sample}/allgenecalled.{sample}_1.coverage.tsv",
        unpaired = "{project}/genecatalog/{sample}/allgenecalled.{sample}_unpaired.coverage.tsv"
    threads: 16
    run:
        shell("source /data/tools/samtools/1.3/env.sh; source /data/tools/BamM/1.7.0/env.sh; bamm parse -c {output.paired} -m counts -b {input.paired} -t {threads}")
        shell("source /data/tools/samtools/1.3/env.sh; source /data/tools/BamM/1.7.0/env.sh; bamm parse -c {output.unpaired} -m counts -b {input.unpaired} -t {threads}")

rule fixcounts:
    input:
        paired = "{project}/genecatalog/{sample}/allgenecalled.{sample}_1.coverage.tsv",
        unpaired = "{project}/genecatalog/{sample}/allgenecalled.{sample}_unpaired.coverage.tsv"  
    output:
        "{project}/genecatalog/{sample}/allgenecalled.{sample}.coverage.tsv",
    run:
        import pandas as pd
        import numpy as np

        t1 = pd.read_table(input.paired)
        t2 = pd.read_table(input.unpaired)

        # Extract the sample name from the filepath of the bam file in the third column
        sample = t1.columns[2].split('.')[-2]

        # The first two columns (gene name and gene length) are the same for both files
        combined = t1[t1.columns[0:2]]

        # Add a new column with the sample as name
        # Divide the paired counts by two so a fragment with two reads is not counted twice
        # Add the counts of the unpaired data
        combined[sample] = np.floor(t1[t1.columns[2]]/2) + t2[t2.columns[2]]
        combined.to_csv(output[0], sep='\t', index=False)

rule combine_counts:
    input:
        expand("{{project}}/genecatalog/{sample}/allgenecalled.{sample}.coverage.tsv", sample=config["data"])
    output:
        "{project}/genecatalog/all.coverage.tsv"
    run:       
        # Each file has 3 columns: gene,length,count. Get the total. 
        nc = len(input)*3
        # Select only every 3rd col containing the counts
        cols = ",".join(list(map(str,range(3,nc+3,3))))        
        shell("paste -d '\\t' {input} | cut -f 1,2,{cols} > {output}")

rule normalize:
    input:
        "{project}/genecatalog/all.coverage.tsv"
    output:
       "{project}/genecatalog/all.coverage.norm.tpm.tsv"
    run:
        import pandas as pd
        df = pd.read_table(input[0])

        # Select only the colums with counts
        values = df[df.columns[2:]]

        #
        # Normalize by gene length
        #

        # Extract the gene lengths
        length_series = df[df.columns[1]]

        # Create a dataframe with the same size as the values
        # Multiple the length column multiple times
        length_df = pd.concat([length_series for n in range(0,len(values.columns))], axis=1)

        # Make sure the columns have the same names
        length_df.columns = values.columns

        # Divide the counts by the gene length
        norm_length = values.divide(length_df)

        #
        # Normalize by depth
        #

        # Sum each column to get the depth and format it as a single row
        depth = pd.DataFrame(norm_length.sum()).transpose()

        # Concatenate all the rows
        depth_df = pd.concat([depth for n in range(0,len(values))], axis=0)

        # Adjust the index values
        depth_df.index = values.index

        # Divide by sequencing depth and multiply with scaling factor
        norm_length_depth = norm_length.divide(depth_df) * 1000000

        # Combine the first two columns of the original data file with the normalized counts 
        final = pd.concat([df[df.columns[:2]],norm_length_depth], axis=1)

        # Write final output to disk
        final.to_csv(output[0], sep='\t', index=False)

# Calculate Average Genome Size
# Can be used for normalization
rule ags_microbecensus:
    input:
        forward="{project}/trimming/{sample}_1.fastq.gz",
        reverse="{project}/trimming/{sample}_2.fastq.gz",
        unpaired="{project}/trimming/{sample}_unpaired.fastq.gz"
    output:
        "{project}/microbecensus/{sample}.ags.txt"
    threads: 16
    run:
        reads = ",".join(input)    
        shell("run_microbe_census.py -t {threads} -n 100000000 -l 300 {reads} {output}")

rule to_hdf5:
    input:
        "{project}/genecatalog/all.coverage.norm.tpm.tsv"
    output:
        "{project}/genecatalog/all.coverage.norm.tpm.biom"
    run: 
       shell("set +u; source /data/tools/biom-format/2.1.5/env.sh; set -u; biom convert -i {input} -o {output} --to-hdf5")

rule diamond_genes:
    input:
        fasta="{project}/genecatalog/allgenecalled.fna.gz"
    output:
        tsv="{project}/genecatalog/allgenecalled.diamond.nr.tsv",
        taxonomy="{project}/genecatalog/allgenecalled.diamond.nr-taxonomy.tsv"
    params:
        reference=config["diamond_database"],
        version="0.7.10",
        output="{project}/genecatalog/allgenecalled.diamond.nr",
        format="tab",
        tmp="/tmp"        
    threads: 16
    run:
        shell("/data/tools/diamond/{params.version}/bin/diamond blastx --sensitive -c 1 -d {params.reference} -t {params.tmp} -p {threads} -q {input} -a {params.output}")
        shell("/data/tools/diamond/{params.version}/bin/diamond view -f {params.format} -a {params.output}.daa -o {output.tsv}")
        shell("java -Xmx32G -Djava.awt.headless=true -Duser.language=en -Duser.region=US -cp '/data/tools/MEGAN/5.10.0/jars/MEGAN.jar:/data/tools/MEGAN/5.10.0/jars/data.jar' megan.tools.TaxonPathClassifier -i {output.tsv} -f Detect -ms 50 -me 0.01 -tp 50 -g2t /data/db/megan/gi_taxid_prot-4March2015.bin -o {output.taxonomy}")

rule kraken_genes:
    input:
        fasta="{project}/genecatalog/allgenecalled.fna.gz"
    output:
        kraken = "{project}/genecatalog/allgenecalled.kraken",
        taxonomy = "{project}/genecatalog/allgenecalled.kraken.taxonomy",
        report = "{project}/genecatalog/allgenecalled.kraken.report"
    log:
        "{project}/genecatalog/allgenecalled.kraken.log"
    threads: 16
    run:
        shell("/data/tools/kraken/0.10.5-beta/bin/kraken --preload --db /data/db/kraken/ --threads {threads} --fasta-input {input} --gzip-compressed --check-names --output {output.kraken} 2> {log}")
        shell("/data/tools/kraken/0.10.5-beta/bin/kraken-translate --db /data/db/kraken/ --mpa-format {output.kraken} > {output.taxonomy}")
        shell("/data/tools/kraken/0.10.5-beta/bin/kraken-mpa-report --db /data/db/kraken/  {output.kraken} > {output.report}")

rule kraken_filter:
    input:
        taxonomy = "{project}/genecatalog/allgenecalled.kraken.taxonomy",
        biom="{project}/genecatalog/all.coverage.tsv"
    output:
        taxonomy = "{project}/genecatalog/allgenecalled.kraken.filtered.taxonomy",
    run:
        # Add missing genes as unclassified
        # The list needs to be complete in order to merge with the BIOM file
        genes = [line.strip().split()[0] for line in open(input.biom).readlines()][1:]

        genedict = {}
        for line in open(input.taxonomy):
            gene, tax = line.split()
            genedict[gene] = tax

        out = open(output.taxonomy, 'w')
        for gene in genes:
            if gene not in genedict:
                genedict[gene] = 'd__unclassified|p__unclassified|c__unclassified|o__unclassified|f__unclassified|g__unclassified|s__unclassified'
            out.write("%s\t%s\n" % (gene,genedict[gene].replace('|',';')))
        out.close()

rule biom_add_kraken_taxonomy:
    input:
        biom="{project}/genecatalog/all.coverage.norm.tpm.biom",
        taxonomy = "{project}/genecatalog/allgenecalled.kraken.filtered.taxonomy",
    output:
       "{project}/genecatalog/all.coverage.norm.tpm.taxonomy.kraken.biom"
    run:
        shell("set +u; source /data/tools/biom-format/2.1.5/env.sh; set -u; biom add-metadata -i {input.biom} -o {output} --observation-metadata-fp {input.taxonomy} --observation-header gene,taxonomy --sc-separated taxonomy")       

rule diamond_genes_filter:
    input:
        taxonomy="{project}/genecatalog/allgenecalled.diamond.nr-taxonomy.tsv",
        biom="{project}/genecatalog/all.coverage.tsv"
    output:
        taxonomy="{project}/genecatalog/allgenecalled.diamond.nr-taxonomy-filtered.tsv"
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
                    level,value = re.match('^\[(.*)\]\s(.*)', s[i].strip()).groups()                    
                    if int(s[i+1]) >= minscore:
                        taxdict[level] = value
                except:
                    pass
            taxonomy = []
            for level in ["SuperKingdom","Phylum","Class","Order","Family","Genus","Species"]:
                taxonomy.append(taxdict.setdefault(level, "unclassified"))

            taxstring = ";".join(taxonomy)
            genestring = ";".join(s[0:1])
            filteredtax[genestring] = taxstring

        # Add missing genes as unclassified
        # The list needs to be complete in order to merge with the BIOM file
        genes = [line.strip().split()[0] for line in open(input.biom).readlines()][1:]        
        for gene in genes:
            if gene not in filteredtax:
                filteredtax[gene] = 'unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified'        
            out.write("%s\t%s\n" % (gene,filteredtax[gene]))
        out.close()

rule diamond_genes_split:
    input:
        taxonomy="{project}/genecatalog/allgenecalled.diamond.nr-taxonomy-filtered.tsv"
    output:
        taxonomy="{project}/genecatalog/allgenecalled.diamond.nr-taxonomy-filtered-split.tsv"
    shell: "cat {input} | tr ';' '\t' > {output}"


rule biom_add_diamond_taxonomy:
    input:
        biom="{project}/genecatalog/all.coverage.norm.tpm.biom",
        taxonomy="{project}/genecatalog/allgenecalled.diamond.nr-taxonomy-filtered.tsv"
    output:
       "{project}/genecatalog/all.coverage.norm.tpm.taxonomy.biom"
    run:
        shell("set +u; source /data/tools/biom-format/2.1.5/env.sh; set -u; biom add-metadata -i {input.biom} -o {output} --observation-metadata-fp {input.taxonomy} --observation-header gene,taxonomy --sc-separated taxonomy")

rule merge_taxonomy:
    input:
        table="{project}/genecatalog/all.coverage.norm.tpm.tsv",
        taxdiamond="{project}/genecatalog/allgenecalled.diamond.nr-taxonomy-filtered.tsv",
        taxkraken="{project}/genecatalog/allgenecalled.kraken.filtered.taxonomy",
    output:
        "{project}/genecatalog/all.coverage.norm.tpm.taxonomy.tsv"
    run:
        import pandas as pd
        df = pd.read_table(input.table)
        df.index = df[df.columns[0]]

        taxdiamond = pd.read_table(input.taxdiamond,header=None, names=['gene','taxdiamond'])
        taxdiamond.index = taxdiamond.gene

        taxkraken = pd.read_table(input.taxkraken,header=None, names=['gene','taxkraken'])
        taxkraken.index = taxkraken.gene


        table_with_taxonomy = pd.concat([df,taxdiamond.taxdiamond, taxkraken.taxkraken],axis=1)

        # Write final output to disk
        table_with_taxonomy.to_csv(output[0], sep='\t', index=False)


#
# Taxator-tk
#
rule taxator_align:
    input:
        fasta="{project}/genecatalog/allgenecalled.fna.gz"
    output:
        alignment="{project}/genecatalog/taxator-tk/R1.alignments.gz"
    threads: 16
    run:
        #shell("awk '{{print $1}}' < {input.forward} > {output.fasta}") 
        shell("set +u; source /data/tools/taxator-tk/1.3.0e-rc2/env.sh; source /data/tools/last/418/env.sh; set -u; zcat {input} | parallel-fasta -j {threads} 'lastal -f 1 /data/tools/taxator-tk/1.2.1-extended/aligner-index/last/nonredundant-microbial_20121122/dna' | lastmaf2alignments | sort | gzip > {output.alignment}")

rule taxator_predict:
    input:
        fasta="{project}/genecatalog/allgenecalled.fna",
        alignment="{project}/genecatalog/taxator-tk/R1.alignments.gz"
    output:
        "{project}/genecatalog/taxator-tk/R1.gff"
    threads: 16
    shell: "export TAXATORTK_TAXONOMY_NCBI=/data/tools/taxator-tk/1.2.1-extended/refdata/nonredundant-microbial_20121122/ncbi-taxonomy/; zcat {input.alignment} | /data/tools/taxator-tk/taxator-tk/Build-x86_64/taxator -a rpa -q {input.fasta} -f /data/tools/taxator-tk/1.2.1-extended/refdata/nonredundant-microbial_20121122/refdata.fna  -g /data/tools/taxator-tk/1.2.1-extended/refdata/nonredundant-microbial_20121122/mapping.tax -p {threads} > {output}"

rule taxator_bin:
    input:
        "{project}/genecatalog/taxator-tk/R1.gff"
    output:
        "{project}/genecatalog/taxator-tk/R1.tax"
    # TODO: set binner settings
    shell: "sort {input} | /data/tools/taxator-tk/taxator-tk/Build-x86_64/binner > {output}"

rule taxator_annotate:
    input:
        "{project}/genecatalog/taxator-tk/R1.tax"
    output: 
        "{project}/genecatalog/taxator-tk/R1.taxonomy"
    shell: "/data/tools/taxator-tk/taxator-tk/Build-x86_64/taxknife -f 2 -m annotate -a -s path  < {input} > {output}"


rule uproc_genes:
    input:
        "{project}/genecatalog/allgenecalled.faa.gz"
    output:
        kegg="{project}/genecatalog/uproc/allgenecalled.uproc.kegg.txt",
        cog="{project}/genecatalog/uproc/allgenecalled.uproc.cog.txt",
        pfam="{project}/genecatalog/uproc/allgenecalled.uproc.pfam.txt",
    threads: 8
    run:
        shell("uproc-prot --preds -o {output.kegg} /data/db/uproc/kegg_20140317/ /data/db/uproc/model {input}")
        shell("uproc-prot --preds -o {output.cog} /data/db/uproc/cog2014/ /data/db/uproc/model {input}")
        shell("uproc-prot --preds -o {output.pfam} /data/db/uproc/pfam28/ /data/db/uproc/model {input}")

rule uproc_filter_ko:
    input:
        kegg="{project}/genecatalog/uproc/allgenecalled.uproc.kegg.txt",
        biom="{project}/genecatalog/all.coverage.tsv"
    output:
        kegg="{project}/genecatalog/uproc/allgenecalled.uproc.kegg.filtered.txt",
    run:
        kegg_dict = {}
        for line in open(input.kegg):
            id, gene, length, ko, score = line.strip().split(',')
            if gene in kegg_dict:
                if kegg_dict[gene]['score'] < score:
                    kegg_dict[gene] = {'ko' : ko, 'score': score}
            else:
                kegg_dict[gene] = {'ko' : ko, 'score': score}

        ko_descr = {'None' : 'None'}
        for line in open('/data/db/uproc/kegg_20140317/ko_description.txt').readlines():
           ko, descr = line.strip().split('\t',1)
           ko_descr[ko.split(':')[1]] = descr

        out = open(output.kegg, 'w')
        genes = [line.strip().split()[0] for line in open(input.biom).readlines()][1:]        

        for gene in genes:
            if gene not in kegg_dict:
                kegg_dict[gene] = {'ko' : 'None'}
            print(ko_descr.setdefault(kegg_dict[gene]['ko'], "Obsolete"))
            out.write("%s\t%s\t%s\n" % (gene, kegg_dict[gene]['ko'], ko_descr.setdefault(kegg_dict[gene]['ko'], "Obsolete")))
        out.close()

rule uproc_filter_pfam:
    input:
        uproc="{project}/genecatalog/uproc/allgenecalled.uproc.pfam.txt"
    output:
        "{project}/genecatalog/uproc/allgenecalled.uproc.pfam.filtered.txt"
    params:
        pfam="/data/db/uproc/pfam28/Pfam-A.clans.tsv.gz"
    run:
        import gzip
        pfams = {}
        for line in gzip.open(params.pfam):
            pfam, empty1, empty2, abbr, description = line.strip().decode("utf-8").split('\t')
            pfams[pfam] = description

        pfam_dict = {}
        for line in open(input.uproc):
            id, gene, length, pfam, score = line.strip().split(',')
            pfam_dict.setdefault(gene,[]).append(pfam)
        
        out = open(output[0], 'w')
        for gene in pfam_dict:
            pfam_list = []
            for pfam in pfam_dict[gene]:
                pfam_list.append("%s;%s" % (pfam, pfams[pfam]))
            out.write("%s\t%s\n" % (gene, "|".join(pfam_list)))
        out.close()

rule add_ko:
    input:
        table="{project}/genecatalog/all.coverage.norm.tpm.taxonomy.tsv",
        ko="{project}/genecatalog/uproc/allgenecalled.uproc.kegg.filtered.txt"
    output:
        "{project}/genecatalog/all.coverage.norm.tpm.taxonomy.ko.tsv",
    run:
        import pandas as pd
        df = pd.read_table(input.table)
        df.index = df[df.columns[0]]

        ko = pd.read_table(input.ko,header=None,names=['gene','ko_best','ko_best_description'])
        ko.index = ko.gene

        final = pd.concat([df,ko[ko.columns[1:]]],axis=1)

        # Write final output to disk
        final.to_csv(output[0], sep='\t', index=False)

rule add_pfam:
    input:
        table="{project}/genecatalog/all.coverage.norm.tpm.taxonomy.ko.tsv",
        annot="{project}/genecatalog/uproc/allgenecalled.uproc.pfam.filtered.txt"
    output:
        "{project}/genecatalog/all.coverage.norm.tpm.taxonomy.ko.pfam.tsv",
    run:
        import pandas as pd
        df = pd.read_table(input.table)
        df.index = df[df.columns[0]]

        annot = pd.read_table(input.annot,header=None,names=['gene','pfam'])
        annot.index = annot.gene

        final = pd.concat([df,annot[annot.columns[1:]]],axis=1)

        # Write final output to disk
        final.to_csv(output[0], sep='\t', index=False)


rule add_taxlevels:
    input:
        table="{project}/genecatalog/all.coverage.norm.tpm.taxonomy.ko.pfam.tsv",
        taxonomy="{project}/genecatalog/allgenecalled.diamond.nr-taxonomy-filtered-split.tsv"
    output:
        table="{project}/genecatalog/all.coverage.norm.tpm.taxonomy.ko.pfam.taxlevels.tsv",
    run:
        import pandas as pd
        df = pd.read_table(input.table)
        df.index = df[df.columns[0]]

        taxdiamond = pd.read_table(input.taxonomy,header=None, names=["gene", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"])
        taxdiamond.index = taxdiamond.gene


        table_with_taxonomy = pd.concat([df,taxdiamond],axis=1)

        # Write final output to disk
        table_with_taxonomy.to_csv(output[0], sep='\t', index=False)

    
rule aggregate_taxonomy_and_ko:
    input:
        table="{project}/genecatalog/all.coverage.norm.tpm.taxonomy.ko.pfam.taxlevels.tsv",
    output:
        table="{project}/genecatalog/all.coverage.norm.tpm.taxonomy.ko.pfam.taxlevels.aggregated.tsv",
    run:
        import pandas as pd
        import numpy as np
        features = pd.read_table(input.table)
        ko_groupby = features.groupby(['ko_best','Family'])
        ko_mean = ko_groupby.aggregate(np.sum)
        ko_mean.to_csv(open(output.table, 'w'),sep='\t')

