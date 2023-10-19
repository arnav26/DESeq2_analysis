README for Kallisto to DESeq2 Processing Pipeline. Written during my PhD in 2019.
# Tximport

This R script provides a pipeline for importing Kallisto quantification results, processing them into gene-level counts, and saving these counts for downstream differential expression analysis (e.g., with DESeq2).
Dependencies:

Ensure that the following R libraries are installed:

    tximport
    biomaRt

Workflow Steps:

    Get Transcript to Gene Mapping: Fetches the mapping between transcripts and genes using the biomaRt package for the mouse genome (Mus musculus) from the Ensembl database.
    Retrieve Kallisto Output Files: Detects the Kallisto output files (abundance.h5) from the working directory.
    Import and Process Kallisto Files: Imports the Kallisto outputs and aggregates the transcript-level quantifications to gene-level counts using the tximport package.
    Save Counts to File: The processed gene counts are saved to a file named gene_counts.

Usage:

    Working Directory Setup:
        Ensure that the R working directory contains folders for each sample, and within each folder, there should be an abundance.h5 file which is the Kallisto output.
        If your directory structure is different, modify the get_kallisto_files function accordingly.

    Output File Produced: gene_counts

    To Run the Script:

    Simply source the script in R, and the pipeline will execute the steps sequentially:

    R

    source("path_to_script.R")

Customization:

    Adjust the biomaRt parameters if you are working with a different organism.
    If your Kallisto output files are in a different directory or have a different naming convention, adjust the get_kallisto_files function.


# DESeq2_analysis
Differential gene expression analysis using DESeq2 using gene counts from Kallisto quantification
README for DESeq2 Analysis Pipeline

This R script provides a pipeline to process and analyze gene expression count data using the DESeq2 package. The pipeline encompasses reading the data, differential expression analysis, plotting PCA and heatmaps, and creating ranked lists for the detected genes based on the p-values and fold changes.
Dependencies:

Ensure that the following R libraries are installed:

    DESeq2
    tidyverse
    vsn
    pheatmap
    RColorBrewer

Workflow Steps:

    Load Libraries: Initializes necessary libraries for the subsequent steps.
    Read and Prepare Data: Reads in the sample data and associated meta-information. The function also filters the data to keep only the genes with counts greater than or equal to 10 in at least 3 samples.
    Variance Transformation: Applies variance stabilizing transformation (VST) to the data, which is useful for downstream applications like PCA and heatmap plotting.
    Plot Heatmap and PCA: Generates a sample-to-sample distance heatmap and a PCA plot which can help in visualizing overall patterns and sample clustering.
    Differential Expression Analysis: Performs differential expression analysis using a given design model and merges with ortholog data.
    Rank Lists and Writing to File: Creates ranked lists based on p-values and fold changes for the detected genes and writes them to the specified files.

Usage:

    Input Files Required:
        newmodel: A table that contains the sample model (meta-information).
        tximport_raw_gene_counts: Raw gene counts table.
        manual_orthologs.csv: Manual ortholog data for the genes.

    Output Files Produced:
        GhurEuroPrerank.rnk
        GhurEuroPostrank.rnk
        GhurPostPrerank.rnk
        EuroPostPrerank.rnk
