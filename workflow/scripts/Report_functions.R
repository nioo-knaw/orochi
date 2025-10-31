
# READ BASED --------------------------------------------------------------

# Loading required packages --------------------------------------------------------
#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager", repos = "https://cloud.r-project.org")

#BiocManager::install("MicrobiotaProcess")

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

# report2 options
metaphlan_secondary <- params$metaphlan_secondary
config_file <- params$config_file

#config_file <- here("config", "configfile.yaml") # change to the relative corresponding path
config <- yaml::read_yaml(config_file)

# merged_table.txt (output from rule MetaPhlAn_secondary)
#metaphlan_secondary <- snakemake@input[["metaphlan_secondary"]]
# samples.tsv FROM CONFIG
samples <- config$samples
# outdir FROM CONFIG
outdir <- config$outdir
# dir where \\.emapper\\.annotations\\.adjusted\\.tsv is (line 423)
# ...

plotsdir <- "results/08_plots/PLOTS/1-Reads"
plotsdirbins <- "results/08_plots/PLOTS/3-Bins"
dir.create(file.path(outdir, plotsdir), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, plotsdirbins), recursive = TRUE, showWarnings = FALSE)

# To save the other htmls in one same location for easier visualization later
dir.create(file.path(outdir, "results/08_plots/rsc"), recursive = TRUE, showWarnings = FALSE)
rsc_path <- file.path(outdir, "results/08_plots/rsc")

# Making object -----------------------------------------------------------
#(https://github.com/YuLab-SMU/MicrobiotaProcess/issues/58)
#mpse <- mp_import_metaphlan(profile='./merged_table.txt') ## CHANGE output from rule MetaPhlAn_secondary

## TO CONNECT TO SNAKEMAKE
# Access input file
mpse <- mp_import_metaphlan(profile = metaphlan_secondary)

# Including pools information
#sample_groups <- read.table("../orochiNIOO/config/inUse_samples.tsv", header = TRUE) ## CHANGE
#colnames(sample_groups)[1] <- "Sample"
#sample_groups <- sample_groups[,c(1,length(sample_groups))]
#mpse3 <- mpse %>%
#  left_join(sample_groups, by = "Sample")

# Including TREATMENT information (in my case, treatment1)
sample_groups <- read.table(samples, header = TRUE)
colnames(sample_groups)[1] <- "Sample"
#colnames(sample_groups)[2] <- "Treatment1" ### Make it so they can choose which treatment column to use?
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

## MAKING ONE BAR PER POOL
# onebar_taxa_plots <- list()
# 
# for (level in taxa_levels) {
#   # Collapse to mean abundance per Family per treatment
#   collapsed_family <- mpse3 %>%
#     group_by(!!sym(level), treatment1, Sample) %>%
#     summarise(taxon_abund = sum(RelAbundanceBySample), .groups = "drop") %>%
#     group_by(!!sym(level), treatment1) %>%
#     summarise(mean_abund = mean(taxon_abund), .groups = "drop")
# 
#   # Reshape: wide format (rows = Family, cols = treatment, values = abundance)
#   mat <- collapsed_family %>%
#     pivot_wider(names_from = treatment1, values_from = mean_abund, values_fill = 0) %>%
#     as.data.frame()
#   rownames(mat) <- mat[,1]
#   mat[,1] <- NULL
# 
#   # Step 2: build MPSE object correctly
#   mpse_treatment <- MPSE(
#     assays = list(Abundance = as.matrix(mat))
#   )
# 
#   # Step 3: plot with mp_plot_abundance
#   p3 <- mpse_treatment %>%
#     mp_plot_abundance(
#       .abundance = Abundance,
#       taxa.class = level,
#       relative = FALSE,
#       force = TRUE,
#       topn = 20
#     )
#   p3 <- p3 + labs(title = paste("(Bacterial) Taxonomic Abundance at", level)) +
#     theme(plot.title = element_text(hjust = 0.5))
#   
#   # Save each plot with its level in filename
#   ggplot2::ggsave(
#     filename = paste0("plots/2-taxonomic_abundance_average_", tolower(level), ".tiff"),
#     plot = p3,
#     width = 12,
#     height = 10,
#     units = "in",
#     dpi = 500,
#     compression = "lzw"
#   )
#   onebar_taxa_plots[[level]] <- p3
# }
# 
# plot2_average_taxonomic_abundance <- function(level = "Family") {
#   if (!level %in% names(onebar_taxa_plots)) {
#     stop("Invalid level! Choose from: ", paste(names(onebar_taxa_plots), collapse = ", "))
#   }
#   onebar_taxa_plots[[level]]
# }

## As heatmaps

taxa_plots_heatmap <- list()
library(ggplotify)

# Loop over each taxonomic level
for (level in taxa_levels) {
  
  h1 <- mpse3 %>%
    mp_plot_abundance(
      .abundance = Abundance,
      #.group = sample_pool,
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
      #.group = sample_pool,
      #.group = treatment1,
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

# The size of point also can be mapped to other variables such as Observe, or Shannon 
# Then the alpha diversity and beta diversity will be displayed simultaneously.
pcoa2 <- mpse3 %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = treatment1, 
    .color = treatment1, 
    .size = 4, #Observe, 
    .alpha = 0.8, #Shannon,
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

contig_files <- list.files(base_dir, pattern = "_linkages_by_contig\\.txt$", 
                           full.names = TRUE, recursive = TRUE)

all_plots <- list()

for (i in seq_along(contig_files)) {
  contig_file <- contig_files[i]
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
  all_plots[[i]] <- plot
}

plot7_maglinkage <- function() {
  for(i in seq_along(all_plots)) {
    print(all_plots[[i]])
  }
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
#NIOO <- paste0(getwd(), "/NIOO.gif")
WUR <- paste0(getwd(), "/WUR.png")

NIOO <- "NIOO.gif"


