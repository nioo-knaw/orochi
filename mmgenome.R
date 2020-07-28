library(mmgenome)
load(snakemake@input$mmgenome)

p <- mmplot(data = d, x = "plate1", y = "plate2", log.x = T, log.y = T, color = "essential", minlength = 3000)
