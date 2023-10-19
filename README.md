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
