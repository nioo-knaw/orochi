library(rmarkdown)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript render_report.R <config_path> <metaphlan_secondary> <output_file>", call. = FALSE)
}

# get paths from snakemake
#rmd_file <- file.path(snakemake@scriptdir, "OROCHIPlots.Rmd")

config_path <- normalizePath(args[1])
metaphlan_secondary <- normalizePath(args[2])
output_path <- normalizePath(args[3], mustWork = FALSE)

# Get script directory
args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_full[grep("--file=", args_full)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()

rmd_file <- file.path(script_dir, "OROCHIPlots.Rmd")

render(
  input = rmd_file,
  output_file = output_path,
  params = list(
    metaphlan_secondary = metaphlan_secondary,
    config_file = config_path
  )
)
