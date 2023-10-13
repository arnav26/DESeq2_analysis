# DESeq2_analysis
Differential gene expression analysis using DESeq2 using gene counts from Kallisto quantification
Step 1: Run tximport.R in the kallisto output directory to compile the gene level aggregate counts.
Step 2. Run DEseq.R with the results of step1 + experiment setup to generate DE results for each comparison.
