library(tximport)
library(biomaRt)

# Function to get transcript to gene mapping from biomart
get_transcript_to_gene <- function() {
  mart <- biomaRt::useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "mmusculus_gene_ensembl",
    host = 'ensembl.org'
  )
  
  t2g <- biomaRt::getBM(
    attributes = c("ensembl_transcript_id_version", "ensembl_gene_id"),
    mart = mart
  )
  
  return(t2g)
}

# Function to get file paths of Kallisto output files
get_kallisto_files <- function() {
  samples <- list.files()
  files <- file.path(samples, "abundance.h5")
  
  return(files)
}

# Function to import and process Kallisto files
import_and_process_kallisto <- function(files, t2g) {
  txi.kallisto <- tximport(
    files,
    type = "kallisto",
    tx2gene = t2g,
    countsFromAbundance = "lengthScaledTPM"
  )
  
  colnames(txi.kallisto$counts) <- basename(files)
  
  return(txi.kallisto)
}

# Function to save counts to file
save_counts_to_file <- function(txi_kallisto) {
  counts <- txi_kallisto$counts
  counts <- round(counts)
  write.table(counts, "gene_counts", quote = FALSE)
}

# Main function to execute steps
run_analysis <- function() {
  t2g <- get_transcript_to_gene()
  kallisto_files <- get_kallisto_files()
  txi_kallisto <- import_and_process_kallisto(kallisto_files, t2g)
  save_counts_to_file(txi_kallisto)
}

# Execute the analysis
run_analysis()
