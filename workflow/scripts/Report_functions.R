
# READ BASED --------------------------------------------------------------

# Loading required packages --------------------------------------------------------

library(MicrobiotaProcess)
library(ggplot2)
library(vegan)
library(ggside)
library(htmltools)
library(dplyr)
library(tidyr)
library(ggplotify)
library(ggtree) #(Not available for this version of R)

library(yaml)

# report options
metaphlan_secondary <- params$metaphlan_secondary
config_file <- params$config_file

config <- yaml::read_yaml(config_file)

# samples.tsv FROM CONFIG
samples <- config$samples

# outdir FROM CONFIG
outdir <- config$outdir


plotsdir <- "results/09_plots/PLOTS/1-Reads"
plotsdirbins <- "results/09_plots/PLOTS/3-Bins"
dir.create(file.path(outdir, plotsdir), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, plotsdirbins), recursive = TRUE, showWarnings = FALSE)

# To save the other htmls in one same location for easier visualization later
#dir.create(file.path(outdir, "results/08_plots/rsc"), recursive = TRUE, showWarnings = FALSE)
rsc_path <- file.path(outdir, "results/09_plots/rsc")

# Making object -----------------------------------------------------------
#(https://github.com/YuLab-SMU/MicrobiotaProcess/issues/58)

## TO CONNECT TO SNAKEMAKE
# Access input file
mpse <- mp_import_metaphlan(profile = metaphlan_secondary)

# Including TREATMENT information (in my case, treatment1)
sample_groups <- read.table(samples, header = TRUE)
colnames(sample_groups)[1] <- "Sample"
sample_groups <- sample_groups[,c(1,2)]
mpse3 <- mpse %>%
  left_join(sample_groups, by = "Sample")

# Alpha diversity analysis ------------------------------------------------

mpse3 %<>%
  mp_cal_abundance( # for each sample
    .abundance = Abundance,
    force = TRUE
  ) %>%
  mp_cal_abundance( # for each group 
    .abundance=Abundance,
    #.group=sample_pool,
    .group=treatment1,
    force = TRUE
  )

# Filter out samples with missing or invalid treatment1 before plotting ## MAKE IT AS FOR treatment2 AS WELL
mpse3 <- mpse3 %>%
filter(!is.na(treatment1) & treatment1 != "")

# Observed and Shannon (per sample)
mpse3 %<>% mp_cal_alpha(.abundance=Abundance, force=TRUE) 
f1 <- mpse3 %>% mp_plot_alpha(.alpha=c(Observe, Shannon))
f1 <- f1 + theme_bw() + labs(title = "Alpha Diversity by Sample") + theme(plot.title = element_text(hjust = 0.5))

# Observed and Shannon (per group)
f2 <- mpse3 %>% mp_plot_alpha(.group=treatment1, .alpha=c(Observe, Shannon))
f2 <- f2 + theme_bw() + labs(title = "Alpha Diversity by Treatment") + theme(plot.title = element_text(hjust = 0.5))
f3 <- f1 / f2

#ggplot2::ggsave(filename = "plots/1-alpha_diversity.tiff", plot = f3, dpi = 500, width = 12, height = 10, units = "in", compression = "lzw")
ggplot2::ggsave(filename = file.path(outdir, plotsdir, "1-alpha_diversity.tiff"), plot = f3, dpi = 500, width = 12, height = 10, units = "in", compression = "lzw")

plot1_alpha_diversity <- function() {
  f3
}
# Taxonomy abundance ------------------------------------------------------

## As bars
# Per group (pool)
# Define taxonomic levels
taxa_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa_plots <- list()

# Loop over each taxonomic level
for (level in taxa_levels) {
  
  p2 <- mpse3 %>%
    mp_plot_abundance(
      .abundance = Abundance,
      #.group = sample_pool,
      .group = treatment1,
      taxa.class = !!sym(level),  # dynamically use level
      topn = 20,
      relative = FALSE,
      force = TRUE
    )
  
  p2 <- p2 + labs(title = paste("(Bacterial) Taxonomic Abundance at", level)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save each plot with its level in filename
  ggplot2::ggsave(
    filename = file.path(outdir,plotsdir,paste0("2-taxonomic_abundance_", tolower(level), ".tiff")),
    plot = p2,
    width = 12,
    height = 10,
    units = "in",
    dpi = 500,
    compression = "lzw"
  )
  taxa_plots[[level]] <- p2
}

plot2_taxonomic_abundance <- function(level = "Family") {
  if (!level %in% names(taxa_plots)) {
    stop("Invalid level! Choose from: ", paste(names(taxa_plots), collapse = ", "))
  }
  taxa_plots[[level]]
}

## As heatmaps

taxa_plots_heatmap <- list()
library(ggplotify)

# Loop over each taxonomic level
for (level in taxa_levels) {
  
  h1 <- mpse3 %>%
    mp_plot_abundance(
      .abundance = Abundance,
      .group = treatment1,
      taxa.class = !!sym(level),
      relative = TRUE,
      topn = 20,
      geom = 'heatmap',
      sample.dist = 'bray',
      sample.hclust = 'average',
      force=TRUE
    )
  
  h1_gg <- as.ggplot(h1)
  h1_gg <- h1_gg + labs(title = paste("(Bacterial) Taxonomic Abundance at", level)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save each plot with its level in filename
  ggplot2::ggsave(
    filename = file.path(outdir, plotsdir, paste0("3-taxonomic_abundance_heatmap_", tolower(level), ".tiff")),
    plot = h1_gg,
    width = 12,
    height = 10,
    units = "in",
    dpi = 500,
    compression = "lzw"
  )
  taxa_plots_heatmap[[level]] <- h1_gg
}

plot3_taxonomic_abundance_heatmap <- function(level = "Family") {
  if (!level %in% names(taxa_plots_heatmap)) {
    stop("Invalid level! Choose from: ", paste(names(taxa_plots_heatmap), collapse = ", "))
  }
  taxa_plots_heatmap[[level]]
}



## MAKING ONE BAR PER POOL
onebar_taxa_plots <- list()
average_taxa_plots_heatmap <- list()

for (level in taxa_levels) {
  # Collapse to mean abundance per Family per treatment
  collapsed_family <- mpse3 %>%
    group_by(!!sym(level), treatment1, Sample) %>%
    summarise(taxon_abund = sum(RelAbundanceBySample), .groups = "drop") %>%
    group_by(!!sym(level), treatment1) %>%
    summarise(mean_abund = mean(taxon_abund), .groups = "drop")
  
  # Reshape: wide format (rows = Family, cols = treatment, values = abundance)
  mat <- collapsed_family %>%
    pivot_wider(names_from = treatment1, values_from = mean_abund, values_fill = 0) %>%
    as.data.frame()
  rownames(mat) <- mat[,1]
  mat[,1] <- NULL
  
  # Step 2: build MPSE object correctly
  mpse_treatment <- MPSE(
    assays = list(Abundance = as.matrix(mat))
  )
  
  # Step 3: plot with mp_plot_abundance (barplot)
  p3 <- mpse_treatment %>%
    mp_plot_abundance(
      .abundance = Abundance,
      taxa.class = level,
      relative = FALSE,
      force = TRUE,
      topn = 20
    )
  p3 <- p3 + labs(title = paste("(Bacterial) Taxonomic Abundance at", level)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save each plot with its level in filename
  ggplot2::ggsave(
    filename = file.path(outdir, plotsdir,paste0("2-taxonomic_abundance_average_", tolower(level), ".tiff")),
    plot = p3,
    width = 12,
    height = 10,
    units = "in",
    dpi = 500,
    compression = "lzw"
  )
  onebar_taxa_plots[[level]] <- p3
  
  # Step 3: plot with mp_plot_abundance (heatmap)
  h2 <- mpse_treatment %>%
    mp_plot_abundance(
      .abundance = Abundance,
      taxa.class = !!sym(level),
      relative = TRUE,
      topn = 20,
      geom = 'heatmap',
      sample.dist = 'bray',
      sample.hclust = 'average',
      force=TRUE
    )
  h2_gg <- as.ggplot(h2)
  h2_gg <- h2_gg + labs(title = paste("(Bacterial) Taxonomic Abundance at", level)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save each plot with its level in filename
  ggplot2::ggsave(
    filename = file.path(outdir, plotsdir, paste0("3-taxonomic_abundance_heatmap_average_", tolower(level), ".tiff")),
    plot = h2_gg,
    width = 12,
    height = 10,
    units = "in",
    dpi = 500,
    compression = "lzw"
  )
  average_taxa_plots_heatmap[[level]] <- h2_gg
  
}

plot2_average_taxonomic_abundance <- function(level = "Family") {
  if (!level %in% names(taxa_plots)) {
    stop("Invalid level! Choose from: ", paste(names(onebar_taxa_plots), collapse = ", "))
  }
  onebar_taxa_plots[[level]]
  
}

plot3_average_taxonomic_abundance_heatmap <- function(level = "Family") {
  if (!level %in% names(average_taxa_plots_heatmap)) {
    stop("Invalid level! Choose from: ", paste(names(average_taxa_plots_heatmap), collapse = ", "))
  }
  average_taxa_plots_heatmap[[level]]
}





# Beta diversity analysis -------------------------------------------------

# Standardization (hellinger comes from here)
mpse3 %<>% 
  mp_decostand(.abundance=Abundance)

## Significance between pools
b3 <- mpse3 %>% mp_cal_dist(.abundance=hellinger, distmethod="bray") %>% mp_plot_dist(.distmethod = bray, .group = treatment1, group.test=TRUE, textsize=2)

ggplot2::ggsave(filename = file.path(outdir, plotsdir, "4-significance_between_pools.tiff"), plot = b3, width = 6, height = 5, units = "in", dpi = 500, compression = "lzw")

plot4_significance_between_pools <- function() {
  b3
}

# PCoA analysis -----------------------------------------------------------

mpse3 %<>% mp_cal_pcoa(.abundance=hellinger, distmethod="bray", .dim = 2)

# We also can perform adonis or anosim to check whether it is significant to the dissimilarities of groups.
#mpse3 %<>% mp_adonis(.abundance=hellinger, .formula=~group, distmethod="bray", permutations=9999, action="add") 
#mpse3 %>% mp_extract_internal_attr(name=adonis)

#pcoa1 <- mpse3 %>%
#  mp_plot_ord(
#    .ord = pcoa, 
#    .group = group, 
#    .color = group, 
#    .size = 1.2,
#    .alpha = 1,
#    ellipse = TRUE,
#    show.legend = FALSE # don't display the legend of stat_ellipse
#  )

pcoa2 <- mpse3 %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = treatment1, 
    .color = treatment1, 
    .size = 4,
    .alpha = 0.8,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse 
  ) +
  scale_size_continuous(
    range=c(0.5, 3),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  )

ggplot2::ggsave(filename = file.path(outdir, plotsdir, "5-PCoA.tiff"), plot = pcoa2, width = 12, height = 10, units = "in", dpi = 500, compression = "lzw")

plot5_pcoa <- function() {
  pcoa2
}

# Hierarchical cluster analysis -------------------------------------------

mpse3 %<>%
  mp_cal_clust(
    .abundance = hellinger, 
    distmethod = "bray",
    hclustmethod = "average",
    action = "add"
  )

sample.clust <- mpse3 %>% mp_extract_internal_attr(name='SampleClust')

p <- ggtree(sample.clust) + 
  geom_tippoint(aes(color=treatment1), size = 8) +
  geom_tiplab(as_ylab = TRUE, size = 14) +
  ggplot2::scale_x_continuous(expand=c(0, 0.01))

ggplot2::ggsave(filename = file.path(outdir, plotsdir, "6-clustering.tiff"), plot = p, width = 12, height = 10, units = "in", dpi = 500, compression = "lzw")

plot6_clustering <- function() {
  p
}


# BINS --------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(ggalluvial)
library(ggnewscale)

base_dir <- file.path(outdir, "results/07_maglinkage")

# Summary of bins per assembly
bin_summary_files <- list.files(file.path(outdir, "results/06_binning/checkm2"), pattern = "quality_report\\.tsv$", 
                                full.names = TRUE, recursive = TRUE)
treatments <- basename(dirname(bin_summary_files))
n_bins <- sapply(bin_summary_files, function(file) length(readLines(file)) - 1)

bin_counts <- data.frame(
  treatment = treatments,
  n_bins = n_bins,
  stringsAsFactors = FALSE
)

# Compute cumulative positions for stacked bar
bin_counts <- bin_counts %>%
  mutate(
    end = cumsum(n_bins),
    start = lag(end, default = 0),
    mid = (start + end) / 2
  )

# Create the plot
bin <- ggplot(bin_counts) +
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 6, fill = treatment)) +
  geom_text(aes(x = mid, y = 3, label = n_bins), color = "black", size = 5, fontface = "bold") +
  scale_fill_brewer(palette = "Oranges") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_fixed(ratio = 0.1) +
  labs(
    x = NULL, y = NULL,
    title = "Total number of bins per treatment (assembly)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )
#ggsave(filename = file.path(outdir, plotsdirbins, paste0(assembly,"_MAGlinkage.tiff")), plot = plot, dpi = 500, width = 12, height = 10, units = "in", compression = "lzw")

plot8_binsummary <- function() {
  print(bin)
}


# MAG-Linkage
contig_files <- list.files(base_dir, pattern = "_linkages_by_contig\\.txt$", 
                           full.names = TRUE, recursive = TRUE)

all_plots <- list()

for (contig_file in contig_files) {
  assembly <- stringr::str_match(contig_file, ".*/([A-Za-z0-9_-]+)/markermag/")[,2]
  genome_file <- file.path(dirname(contig_file), paste0(assembly, "_linkages_by_genome.txt"))
  if (!file.exists(genome_file)) next
  
  message("Processing assembly: ", assembly)
  
  # --- Load data ---
  df1 <- read.table(contig_file, header = TRUE, sep = "\t", quote = "", comment.char = "")
  df2 <- read.table(genome_file, header = TRUE, sep = "\t", quote = "", comment.char = "")
  
  # --- Prepare data ---
  df1 <- df1 %>%
    tidyr::separate(Marker___Genome.total., into = c("MarkerGene", "GenomicSeq_total"), sep = "___") %>%
    mutate(GenomicSeq = sub("\\(.*\\)", "", GenomicSeq_total))
  
  merged <- df1 %>%
    left_join(df2, by = c("MarkerGene", "GenomicSeq"))
  
  merged <- merged %>%
    mutate(
      MarkerGene = factor(MarkerGene, levels = unique(MarkerGene)),
      Contig = factor(Contig, levels = unique(Contig)),
      GenomicSeq = factor(GenomicSeq, levels = unique(GenomicSeq))
    )
  
  plot <- ggplot(merged, aes(axis1 = MarkerGene, axis2 = Contig, axis3 = GenomicSeq, y = Linkage)) +
    geom_alluvium(aes(fill = MarkerGene), width = 1/12, show.legend = FALSE) +
    scale_fill_brewer(palette="Blues") +
    ggnewscale::new_scale_fill() +
    geom_stratum(aes(fill = Round), width = 0.1) +
    scale_fill_manual(values = c("Rd1" = "grey90", "Rd2" = "grey50")) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, nudge_x = 0.08, hjust = 0) +
    scale_x_discrete(limits = c("MarkerGene", "Contig", "GenomicSeq"),
                     labels = c("Marker", "Contig", "Genome")) +
    theme_classic(base_size = 12) +
    labs(title = paste0("16S-MAG Linkage for pool ", assembly)) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(), 
      axis.line = element_blank(),
      legend.position.inside = c(1, 1),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "transparent", color = NA)
    )
  ggplot2::ggsave(filename = file.path(outdir, plotsdirbins, paste0(assembly,"_MAGlinkage.tiff")), plot = plot, dpi = 500, width = 12, height = 10, units = "in", compression = "lzw")
  all_plots[[assembly]] <- plot
}

plot7_maglinkage <- function(assembly) {
  if (!assembly %in% names(all_plots)) {
    stop("Invalid assembly! Available: ", paste(names(all_plots), collapse = ", "))
  }
  all_plots[[assembly]]
}

# For the report ----------------------------------------------------------

library(DT)
library(patchwork)
library(htmltools)
library(tidyr)
library(knitr)

safe_call <- function(fun_name, ...) {
  if (exists(fun_name, mode = "function")) {
    do.call(fun_name, list(...))
  } else {
    message(sprintf("Function %s() is not available in this run.", fun_name))
    return(invisible(NULL))
  }
}

logo_file <- paste0(getwd(), "/Orochi_logo.png")
WUR <- paste0(getwd(), "/WUR.png")
NIOO <- "NIOO.gif"


