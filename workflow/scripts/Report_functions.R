
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
plotsdircontigs <- "results/09_plots/PLOTS/2-Contigs"
plotsdirbins <- "results/09_plots/PLOTS/3-Bins"
dir.create(file.path(outdir, plotsdir), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, plotsdircontigs), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, plotsdirbins), recursive = TRUE, showWarnings = FALSE)

# To save the other htmls in one same location for easier visualization later
#dir.create(file.path(outdir, "results/08_plots/rsc"), recursive = TRUE, showWarnings = FALSE)
rsc_path <- file.path(outdir, "results/09_plots/rsc")
dir.create(rsc_path, recursive = TRUE, showWarnings = FALSE)

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

# CONTIGS --------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(httr)

base_dir_05 <- file.path(outdir, "results/05_prokaryote_annotation")
treatments <- list.dirs(file.path(base_dir_05, "eggnog"),
                        full.names = FALSE,
                        recursive = FALSE)

# To get nice names for the KEGG pathways in my data
kegg_map <- read.delim(url("https://rest.kegg.jp/list/pathway/ko"), header = FALSE, stringsAsFactors = FALSE)
colnames(kegg_map) <- c("KEGG_Pathway", "Pathway_name")

# Getting all the pathways to then exclude some
resp <- GET("https://rest.kegg.jp/list/pathway")
txt <- content(resp, "text")
lines <- strsplit(txt, "\n")[[1]]

df <- data.frame(
  pathway_id = str_extract(lines, "map\\d+"),
  pathway_name = str_extract(lines, "(?<=\\t).*"),
  stringsAsFactors = FALSE
)

df <- df %>%
  mutate(category = case_when(
    str_detect(pathway_id, "^map05") ~ "5",
    str_detect(pathway_id, "^map06") ~ "6",
    str_detect(pathway_id, "^map07") ~ "7",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(category))

# To get reference for COG
cog_map <- read.delim(url("https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.fun.tab"), header = FALSE, stringsAsFactors = FALSE)
cog_map <- cog_map[,c(1,4)]
cog_map <- cog_map[cog_map$V4 != "",]
colnames(cog_map) <- c("COG_category", "Description")
write.csv(cog_map, file.path(outdir,plotsdircontigs,"COG_Categories_reference.csv"), row.names = FALSE)
##### THIS ONE CAN BE USED TO PRODUCE THE REFERENCE TABLE THAT WILL BE ADDED IN THE REPORT

merged_list <- list()
kegg_plots_list <- list()
cog_plots_list <- list()
for (trt in treatments) {
  
  salmon_file <- file.path(base_dir_05, "salmon", trt, paste0(trt, "_ORF_TPM.tsv"))
  eggnog_file <- file.path(base_dir_05, "eggnog", trt, paste0(trt, ".emapper.annotations.adjusted"))
  
  # Check files exist (prevents crashing)
  if (file.exists(salmon_file) & file.exists(eggnog_file)) {
    
    salmon_counts <- read.delim(salmon_file, header = TRUE, stringsAsFactors = FALSE)
    eggnog_out   <- read.delim(eggnog_file, header = TRUE, stringsAsFactors = FALSE)
    
    merged_list[[trt]] <- eggnog_out %>%
      inner_join(salmon_counts, by = c("X.query" = "Name"))
    
    kegg_counting <- merged_list[[trt]] %>%
      select(X.query, KEGG_Pathway, all_of(names(merged_list[[trt]])[(which(names(merged_list[[trt]]) == "PFAMs") + 1):ncol(merged_list[[trt]])])) %>%
      filter(KEGG_Pathway != "-") %>%  # remove genes with no KEGG annotation
      mutate(
        KEGG_Pathway = sapply(
          str_split(KEGG_Pathway, ","),
          function(x) paste(x[str_starts(x, "ko")], collapse = ",")
        )
      ) %>%
      filter(KEGG_Pathway != "") %>%
      mutate(KEGG_Pathway = word(KEGG_Pathway, 1, sep = ","))
    
    # Incorporating nice names
    kegg_counting_named <- kegg_counting %>%
      left_join(kegg_map, by = "KEGG_Pathway") %>%
      mutate(Pathway_name = ifelse(is.na(Pathway_name), KEGG_Pathway, Pathway_name)) # fallback if missing
    
    # To get a table with the pathways corresponding to the categories that should be deleted (5-7) as well
    # Removing from my own data the rows that contain a pathway listed in the 5to7 table
    kegg_counting_named_filtered <- kegg_counting_named[!kegg_counting_named$Pathway_name %in% df$pathway_name,]
    
    # Keeping only the further used columns
    kegg_counting_named_adj <- kegg_counting_named_filtered[,3:length(kegg_counting_named_filtered)]
    
    # Average per Pathway for all samples
    kegg_counting_averaged <- kegg_counting_named_adj %>%
      group_by(Pathway_name) %>%
      summarise(across(-last_col(), \(x) mean(x, na.rm = TRUE)))
    
    colnames(kegg_counting_averaged)[1] <- "KEGG Pathway name"
    write.csv(kegg_counting_averaged, file.path(outdir, plotsdircontigs, paste0(trt, "_KEGG_functions.csv")), row.names = FALSE)
    
    ## FOR PLOTTING
    # Keep top 30 for plot
    # First, a column where the sum of all columns is stored to be able to know which ones are the most abundant ones
    kegg_counting_averaged$sum <- rowSums(kegg_counting_averaged[, sapply(kegg_counting_averaged, is.numeric)])
    top_n <- 30
    kegg_count_top <- kegg_counting_averaged %>%
      slice_max(order_by = sum, n = top_n)
    kegg_count_top <- kegg_count_top[,-length(kegg_count_top)]
    
    kegg_count_top_long <- pivot_longer(kegg_count_top, cols = colnames(kegg_count_top[,2:length(kegg_count_top)]), names_to = "Sample", values_to = "Count")
    colnames(kegg_count_top_long)[1] <- "Pathway_name"
    
    # Filtering the ones with a count of zero so they do not show up in the plot
    kegg_count_top_long_filt <- kegg_count_top_long |> 
      dplyr::filter(Count > 0)
    
    # Plot
    k1 <- ggplot(kegg_count_top_long_filt, aes(x = Sample, y = reorder(Pathway_name, Count))) +
      geom_point(aes(size = Count, color = Count)) +
      scale_color_gradient(low = "lightblue", high = "darkblue") +
      theme_minimal() +
      labs(
        x = "",
        y = "KEGG Pathway",
        title = paste0("KEGG Pathway Dotplot (Top 30) - Assembly ", trt),
        color = "TPM",
        size = "TPM"
      )
    
    ggplot2::ggsave(filename = file.path(outdir,plotsdircontigs,paste0(trt,".KEGG_top30.tiff")), plot = k1, width = 8, height = 8, units = "in", dpi = 500, compression = "lzw")
    kegg_plots_list[[trt]] <- k1
    
    
    # Now, COG
    cog_counting <- merged_list[[trt]] %>%
      select(X.query, COG_category, all_of(names(merged_list[[trt]])[(which(names(merged_list[[trt]]) == "PFAMs") + 1):ncol(merged_list[[trt]])])) %>%
      filter(COG_category != "") # remove genes with no COG annotation
    
    # Average per Pathway for all samples
    cog_counting$X.query <- NULL
    cog_counting <- cog_counting[cog_counting$COG_category != "-",]
    cog_counting_averaged <- cog_counting %>%
      group_by(COG_category) %>%
      summarise(across(-last_col(), \(x) mean(x, na.rm = TRUE)))
    
    colnames(cog_counting_averaged)[1] <- "COG Category"
    write.csv(cog_counting_averaged, file.path(outdir, plotsdircontigs, paste0(trt, "_COG_functions.csv")), row.names = FALSE)
    
    ## For the COG plot
    # Keep top 30 for plot
    # First, a column where the sum of all columns is stored to be able to know which ones are the most abundant ones
    cog_counting_averaged$sum <- rowSums(cog_counting_averaged[, sapply(cog_counting_averaged, is.numeric)])
    top_n <- 30
    cog_count_top <- cog_counting_averaged %>%
      slice_max(order_by = sum, n = top_n)
    cog_count_top <- cog_count_top[,-length(cog_count_top)]
    
    cog_count_top_long <- pivot_longer(cog_count_top, cols = colnames(cog_count_top[,2:length(kegg_count_top)]), names_to = "Sample", values_to = "Count")
    colnames(cog_count_top_long)[1] <- "COG_Category"
    
    # Filtering the ones with a count of zero so they do not show up in the plot
    cog_count_top_long_filt <- cog_count_top_long |> 
      dplyr::filter(Count > 0)
    
    # Plot
    c1 <- ggplot(cog_count_top_long_filt, aes(x = Sample, y = reorder(COG_Category, Count))) +
      geom_point(aes(size = Count, color = Count)) +
      scale_color_gradient(low = "#EEE0E5", high = "magenta4") +
      theme_minimal() +
      labs(
        x = "",
        y = "COG Category",
        title = paste0("COG Category Dotplot (Top 30) - Assembly ", trt),
        color = "TPM",
        size = "TPM"
      )
    
    ggplot2::ggsave(filename = file.path(outdir,plotsdircontigs,paste0(trt,".COG_top30.tiff")), plot = c1, width = 6, height = 8, units = "in", dpi = 500, compression = "lzw")
    cog_plots_list[[trt]] <- c1
  }
}


plot1_top30_kegg <- function(trt) {
  kegg_plots_list[[trt]]
}

plot2_top30_cog <- function(trt) {
  cog_plots_list[[trt]]
}

# BINS --------------------------------------------------------------------

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

# Bin scatterplot (it used to be an HTML generated with Python)
# Load libraries
library(readr)
library(plotly)
library(dplyr)
library(tidyr)
library(purrr)
library(htmlwidgets)

bat_base <- file.path(outdir,"/results/06_binning/BAT")
treatments <- list.dirs(bat_base, full.names = FALSE, recursive = FALSE)
base_dir <- file.path(outdir,"/results/06_binning")

all_bin_plots_html <- list()
all_bin_plots_png <- list()
all_bin_plots_html_paths <- list()
all_bin_plots_html_paths_relative <- list()

for (tr in treatments) {
  
  # Construct paths
  bins_path <- file.path(base_dir, "drep/checkm2_genomeinfo", paste0(tr, "_genomeinfo.tsv"))
  tax_path  <- file.path(base_dir, "BAT", tr, paste0(tr, ".bin2classification.names.txt"))
  
  # Load files
  bins <- readr::read_csv(bins_path)
  tax  <- read.delim(tax_path, sep = "\t")
  
  # Example: print names to check
  message("Loaded treatment: ", tr)
  message("  bins: ", bins_path)
  message("  tax:  ", tax_path)
  
  tax <- tax %>%
    mutate(
      taxonomy_concat = pmap_chr(
        select(., superkingdom, phylum, class, order, family, genus, species),
        ~ {
          values <- list(...)
          col_names <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
          
          out <- map2_chr(values, col_names, function(val, col) {
            # Skip "no support" or empty values
            if (is.na(val) || val == "no support") return("")  
            
            prefix <- substr(col, 1, 1)
            paste0(prefix, "_", val)
          })
          
          # Remove empty strings before collapsing
          out <- out[out != ""]
          
          paste(out, collapse = ";")
        }
      )
    )
  binstax <- cbind(bins,tax$taxonomy_concat)
  colnames(binstax)[4] <- "taxonomy_concat"
  
  # Static plot
  p_gg <- ggplot(binstax, aes(x = completeness, y = contamination, color = taxonomy_concat)) +
    geom_point(alpha = 0.7, size = 3) +
    scale_color_brewer(palette = "Set3") +
    theme_minimal() +
    labs(
      title = paste0("Bin Quality (Completeness vs Contamination) - Pool ", tr),
      x = "Completeness (%)",
      y = "Contamination (%)",
      color = "Taxonomy"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      legend.text = element_text(size = 8)
    )
  ggplot2::ggsave(filename = file.path(outdir, plotsdirbins, paste0(tr,"_BinsQuality.tiff")), plot = p_gg, dpi = 500, width = 12, height = 6, units = "in", compression = "lzw")
  all_bin_plots_png[[tr]] <- p_gg
  
  
  # HTML plotly object
  fig <- plot_ly(
    data = binstax,
    x = ~completeness,
    y = ~contamination,
    type = 'scatter',
    mode = 'markers',
    color = ~taxonomy_concat,
    colors = "Set3",
    customdata = ~taxonomy_concat,
    hovertemplate = paste(
      "<b>Completeness:</b> %{x}%<br>",
      "<b>Contamination:</b> %{y}%<br>",
      "<b>Taxonomy:</b> %{customdata}<br>",
      "<extra></extra>"),
    marker = list(
      size = 14,
      opacity = 0.85,
      line = list(width = 0.8, color = 'black')
    )
  ) %>%
    layout(
      title = list(text = paste0("Bin Quality (Completeness vs Contamination) - Pool ", tr), font = list(size=20, family="Arial", color="black")),
      xaxis = list(title = "Completeness (%)", range = c(0,101), showline = TRUE, linecolor = "black", showgrid = TRUE, gridcolor = "lightgrey", zeroline = FALSE),
      yaxis = list(title = "Contamination (%)", range = c(0,100), showline = TRUE, linecolor = "black", showgrid = TRUE, gridcolor = "lightgrey", zeroline = FALSE),
      paper_bgcolor = "white",
      plot_bgcolor = "white",
      font = list(size=14, family="Arial", color="black"),
      legend = list(title = list(text="Taxonomy"))
    )
  html_file <- file.path(outdir, plotsdirbins, paste0(tr,"_BinsQuality.html"))
  saveWidget(fig, html_file, selfcontained = TRUE)
  all_bin_plots_html[[tr]] <- fig
  all_bin_plots_html_paths[[tr]] <- html_file
  # Saving it in another location as well for visualization later
    # Most of the other HTMLs are simply pasted there beforehand, but this one is being created here
  html_file <- file.path(rsc_path, tr, paste0(tr,"_BinsQuality.html"))
  saveWidget(fig, html_file, selfcontained = TRUE)
  all_bin_plots_html_paths_relative[[tr]] <- file.path("rsc", tr, basename(html_file))
}

plot9_binscatterplot <- function(tr) {
  if (!tr %in% names(all_bin_plots_png)) {
    stop("Invalid treatment! Available: ", paste(names(all_bin_plots_png), collapse = ", "))
  }
  all_bin_plots_png[[tr]]
}

plot9_binscatterplot_link <- function(tr) {
  if (!tr %in% names(all_bin_plots_html_paths)) {
    stop("Invalid treatment! Available: ", paste(names(all_bin_plots_html_paths), collapse = ", "))
  }
  html_file <- all_bin_plots_html_paths[[tr]]
  paste0("[View Plotly plot](", html_file, ")")
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


