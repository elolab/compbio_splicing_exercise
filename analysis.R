library(tidyr)
library(pheatmap)

#--------------------------------------------------
# Explore rMATS results on R
#--------------------------------------------------

# What results files there are? 
list.files("rmats_out")


# Read in file related to skipped exon (SE) results using both junction reads and reads covering the targeted exon (JCEC)
rmats.se <- read.table("rmats_out/SE.MATS.JCEC.txt", header=TRUE)
dim(rmats.se)
View(rmats.se)


# Collect PSI values (IncLevel1, IncLevel2) into a matrix
# Separate replicates into individual columns
# e.g. tidyr::separate


# Principal component analysis (PCA) for PSI values
# e.g. prcomp


# Generate a Volcano plot of rMATS results
# x-axis: difference or fold change
# y-axis: -log10 significance


# Identify significantly differential events
# Use the same thresholds as the original study
# (i.e. FDR < 0.05 and absolute deltaPSI > 0.05)


# Order the selected events by absolute difference from largest to smallest.


# Filter the results by coverage
# Split comma-separated read counts to different columns
# For each replicate: inclusion junction count (IJC) + skipping junction count (SJC) > threshold


# Generate a Heatmap of the top 20 splicing events
# e.g. pheatmap::pheatmap


# Plot PSI values for each significant splicing event of Lef1 gene.
# e.g. boxplot (KO, WT)

