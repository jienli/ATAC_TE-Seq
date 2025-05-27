#!/usr/bin/env Rscript

# Load DESeq2
suppressMessages(library("DESeq2"))

# Parse command-line arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3) {
  stop("Usage: Rscript differential_accessibility.R counts.tsv metadata.csv deseq2_results.csv")
}
counts_file <- args[1]
meta_file   <- args[2]
out_file    <- args[3]

# Read count matrix; skip potential comment line from featureCounts
raw_counts <- read.table(counts_file, header=TRUE, sep="\t", comment.char="#", check.names=FALSE)
# Extract only sample count columns (columns 7 onward) and set row names
count_data <- as.matrix(raw_counts[, 7:ncol(raw_counts)])
rownames(count_data) <- raw_counts$Geneid

# Read metadata (sample IDs as row names, must include 'group' column)

col_data <- read.csv(meta_file, header=TRUE, row.names=1, check.names=FALSE)
# Set the first group in metadata as the baseline
col_data$group <- factor(col_data$group)
ref_group <- as.character(col_data$group[1])
col_data$group <- relevel(col_data$group, ref = ref_group)

# Ensure count columns match metadata rows
if(!all(colnames(count_data) %in% rownames(col_data))) {
  stop("Column names of count matrix don't match row names of metadata.")
}
col_data <- col_data[colnames(count_data), , drop=FALSE]

# Create DESeq2 dataset and run differential analysis
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData   = col_data,
                              design    = ~ group)
dds <- DESeq(dds)

# Extract results
res <- results(dds)
res_df <- as.data.frame(res)
res_df$peak <- rownames(res_df)
# Reorder columns
res_df <- res_df[, c("peak", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# Write out CSV
write.csv(res_df, file = out_file, row.names = FALSE)
