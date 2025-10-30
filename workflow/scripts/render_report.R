library(rmarkdown)

# get paths from snakemake
rmd_file <- file.path(snakemake@scriptdir, "OROCHIPlots.Rmd")

config_path <- normalizePath(args[1])
render(
  input = rmd_file,
  output_file = snakemake@output[[1]],
  params = list(
    metaphlan_secondary = snakemake@input[["metaphlan_secondary"]],
    #config_file = snakemake@input[["config"]]
    config_file = config_path
  )
)
