#!/usr/bin/env Rscript

# Script to perform TPM normalization using a table with read counts.
# Usage: Rscript TPM_normalization.R table_raw_counts.txt

# Read count table with the information of all libraries.

args <- commandArgs(trailingOnly = T)

fileName <- args[1]

# Load data

all_reads <- read.table(fileName, header = T)

# Calculate TPM values.

# Step 1. Divide read counts by gene length in kilobases (RPK).

rpk <- all_reads[,3:ncol(all_reads)] / (all_reads$Length/1000)

# Step 2. Sum up RPK values by sample and divide by 1,000,000 (scaling factors).

scalingFactors <- colSums(rpk) / 1000000

# Step 3. Divide each RPK value by its corresponding scaling factor (TPM).

tpm <- sweep(rpk, 2, scalingFactors, FUN = '/')

# Include additional gene information.

tpm <- cbind(all_reads[,1:2], tpm)

# Write annotated TPM values.

write.table(tpm, "TPM_normalized_counts_genes.txt", quote = F, row.names = F, sep = "\t")


