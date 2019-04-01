quast <- read.delim(snakemake@input$quast)
flagstat <- read.delim(snakemake@input$flagstat)
save.image(file=snakemake@output$rdata)

