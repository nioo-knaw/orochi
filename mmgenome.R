library(mmgenome)
ess <- read.table(snakemake@input$essential, header = T, sep = " ")
coverage = read.csv(snakemake@input$coverage, header=T,sep="\t")

var_names <- colnames(coverage)[3:ncol(coverage)]
for(i in var_names) {
  assign(i, data.frame(coverage$X.contig, coverage[[i]]))
}

assembly <- readDNAStringSet(snakemake@input$assembly, format = "fasta")
tax <- read.table(snakemake@input$tax, header = T, sep = "\t")
d <- mmload(assembly = assembly, coverage = c("Shotgun_EPS.bamm.megahit.all.longreads.assembly.plate1_forward_paired.bam"), essential=ess, tax=tax, tax.expand = "Proteobacteria", tax.freq = 85)
rm(list = c("assembly"))
save.image(file=snakemake@output[0])

