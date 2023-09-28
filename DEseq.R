# Load libraries
load_libraries <- function() {
  library(DESeq2)
  library(tidyverse)
  library(vsn)
  library(pheatmap)
  library(RColorBrewer)
}

# Read and prepare data
read_and_prepare_data <- function() {
  model <- read.table("newmodel", header = TRUE)
  levels(model$condition) <- c("Unvaccinated", "Vaccinated")
  counts <- as.matrix(read.table("tximport_raw_gene_counts", header = TRUE, row.names = 'Genes')) %>% round()
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = model,
    design = ~ condition + breed
  )
  
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  
  return(dds)
}

# Variance transformation for PCA and heatmap
variance_transform <- function(dds) {
  vsd <- vst(dds, blind = FALSE)
  return(vsd)
}

# Plot Heatmap and PCA
plot_heatmap_and_pca <- function(vsd) {
  # Heatmap plot
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors
  )
  
  # PCA Plot
  plotPCA(vsd, intgroup = c("breed", "condition"))
}

# DEG Analysis
deg_analysis <- function(dds) {
  dds$group <- as.factor(paste0(dds$condition, dds$breed))
  design(dds) <- ~ group
  dds <- DESeq(dds)
  
  orth <- read.table("manual_orthologs.csv", sep = "\t")
  colnames(orth)[2] <- "Genes"
  
  list(dds = dds, orth = orth)
}

# Rank Lists and Writing to file
rank_lists_and_write <- function(deg_analysis_list) {
  dds <- deg_analysis_list$dds
  orth <- deg_analysis_list$orth
  
  contrasts <- list(
    GhurEuroPre = c("group", "UnvaccinatedGHUR", "UnvaccinatedEURO"),
    GhurEuroPost = c("group", "VaccinatedGHUR", "VaccinatedEURO"),
    GhurPostPre = c("group", "VaccinatedGHUR", "UnvaccinatedGHUR"),
    EuroPostPre = c("group", "VaccinatedEURO", "UnvaccinatedEURO")
  )
  
  rank_files <- c("GhurEuroPrerank.rnk", "GhurEuroPostrank.rnk", "GhurPostPrerank.rnk", "EuroPostPrerank.rnk")
  
  for(i in seq_along(contrasts)) {
    contrast_name <- names(contrasts)[i]
    contrast <- contrasts[[contrast_name]]
    results_df <- cbind(orth$Genes, as.data.frame(results(dds, contrast = contrast)))
    colnames(results_df) <- c("Genes", "log2FoldChange", "padj")
    
    rank_df <- results_df %>% 
      filter(Genes != "") %>% 
      group_by(Genes) %>% 
      mutate(rank = -log10(padj) * sign(log2FoldChange)) %>% 
      select(Genes, rank) %>% 
      filter(!is.na(rank)) %>% 
      arrange(rank)
    
    write.table(rank_df, file = rank_files[i], quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  }
}

# Main function to run the analysis
run_analysis <- function() {
  load_libraries()
  dds <- read_and_prepare_data()
  vsd <- variance_transform(dds)
  plot_heatmap_and_pca(vsd)
  deg_analysis_list <- deg_analysis(dds)
  rank_lists_and_write(deg_analysis_list)
}

# Execute the analysis
run_analysis()
