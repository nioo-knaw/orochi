library(mmgenome)
ess <- read.table(snakemake@input$essential, header = T, sep = " ")


for(i in snakemake@input$coverage) {
  coverage = read.csv(i,header=T,sep="\t")
  assign(colnames(coverage)[3], data.frame(coverage$contig, coverage[,3]))
}

assembly <- readDNAStringSet(snakemake@input$assembly, format = "fasta")
tax <- read.table(snakemake@input$tax, header = T, sep = "\t")

d <- mmload(assembly = assembly, coverage = snakemake@params$samples, essential=ess, tax=tax, tax.expand = "Proteobacteria", tax.freq = 85)
rm(list = c("assembly"))
save.image(file=snakemake@output$mmgenome)

