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
View(rmats.se[,c("IncLevel1","IncLevel2")])
rmats.se <- separate(rmats.se, col="IncLevel1", into=c("PSI_S1R1","PSI_S1R2","PSI_S1R3"), sep=',', remove=TRUE, convert=TRUE)
rmats.se <- separate(rmats.se, col="IncLevel2", into=c("PSI_S2R1","PSI_S2R2","PSI_S2R3"), sep=',', remove=TRUE, convert=TRUE)
View(rmats.se)

# Principal component analysis (PCA) for PSI values
# e.g. prcomp
psi <- rmats.se[,grep("PSI", colnames(rmats.se))]
View(psi)

psi <- na.omit(psi)
psi.pca <- prcomp(t(psi))
plot(psi.pca$x[,"PC1"], psi.pca$x[,"PC2"], pch=21, bg=c(rep("red",3),rep("cyan",3)), xlab="PC1", ylab="PC2", cex=2)
legend("topright", fill=c("red", "cyan"), c("KO", "WT"), bty="n")

summary(psi.pca)$importance
barplot(summary(psi.pca)$importance["Proportion of Variance",], ylab="Proportion of variance explained")

# Generate a Volcano plot of rMATS results
# x-axis: difference or fold change
# y-axis: -log10 significance
plot(rmats.se[,"IncLevelDifference"], -log10(rmats.se[,"FDR"]), pch=21, bg = "gray", xlab="deltaPSI", ylab="-log10FDR", cex=0.5)
abline(h = -log10(0.05), lty=2)
abline(v = 0.05, lty=2)
abline(v = -0.05, lty=2)

# Identify significantly differential events
# Use the same thresholds as the original study
# (i.e. FDR < 0.05 and absolute deltaPSI > 0.05)
rmats.se.filtered <- rmats.se[rmats.se[,"FDR"]<0.05 & abs(rmats.se[,"IncLevelDifference"])>0.05,]
dim(rmats.se.filtered)

# Order the selected events by absolute difference from largest to smallest.
rmats.se.filtered <- rmats.se.filtered[order(abs(rmats.se.filtered[,"IncLevelDifference"]), decreasing = TRUE),]
View(rmats.se.filtered)

# Filter the results by coverage
# Split comma-separated read counts to different columns
# For each replicate: inclusion junction count (IJC) + skipping junction count (SJC) > threshold
rmats.se.filtered <- separate(rmats.se.filtered, col="IJC_SAMPLE_1", into=c("IJC_S1R1","IJC_S1R2","IJC_S1R3"), sep=',', remove=T, convert=T)
rmats.se.filtered <- separate(rmats.se.filtered, col="SJC_SAMPLE_1", into=c("SJC_S1R1","SJC_S1R2","SJC_S1R3"), sep=',', remove=T, convert=T)
rmats.se.filtered <- separate(rmats.se.filtered, col="IJC_SAMPLE_2", into=c("IJC_S2R1","IJC_S2R2","IJC_S2R3"), sep=',', remove=T, convert=T)
rmats.se.filtered <- separate(rmats.se.filtered, col="SJC_SAMPLE_2", into=c("SJC_S2R1","SJC_S2R2","SJC_S2R3"), sep=',', remove=T, convert=T)
View(rmats.se.filtered)

threshold <- 20
sel <- which(rmats.se.filtered[,"IJC_S1R1"] + rmats.se.filtered[,"SJC_S1R1"] > threshold &
               rmats.se.filtered[,"IJC_S1R2"] + rmats.se.filtered[,"SJC_S1R2"] > threshold &
               rmats.se.filtered[,"IJC_S1R3"] + rmats.se.filtered[,"SJC_S1R3"] > threshold &
               rmats.se.filtered[,"IJC_S2R1"] + rmats.se.filtered[,"SJC_S2R1"] > threshold &
               rmats.se.filtered[,"IJC_S2R2"] + rmats.se.filtered[,"SJC_S2R2"] > threshold &
               rmats.se.filtered[,"IJC_S2R3"] + rmats.se.filtered[,"SJC_S2R3"] > threshold)
rmats.se.filtered <- rmats.se.filtered[sel,]
View(rmats.se.filtered)

# Generate a Heatmap of the top 20 splicing events
# e.g. pheatmap::pheatmap
pheatmap(rmats.se.filtered[1:20, grep("PSI",colnames(rmats.se.filtered)),], labels_row=rmats.se.filtered[1:20,"geneSymbol"])

color <- colorRampPalette(c("royalblue","yellow","indianred"))(10)
pheatmap(rmats.se.filtered[1:20, grep("PSI",colnames(rmats.se.filtered)),], labels_row=rmats.se.filtered[1:20,"geneSymbol"], color=color)

# Plot PSI values for each significant splicing event of Lef1 gene.
# e.g. boxplot (KO, WT)
rmats.se.filtered[rmats.se.filtered[,"geneSymbol"]=="Lef1",]

selected.event <- list(KO=unlist(rmats.se.filtered["10187",c("PSI_S1R1","PSI_S1R2","PSI_S1R3")]), WT=unlist(rmats.se.filtered["10187",c("PSI_S2R1","PSI_S2R2","PSI_S2R3")]))
boxplot(selected.event, col=c("red","cyan"), ylab="deltaPSI", main="Lef1", las=1, ylim=c(0,1))
set.seed(1234)
stripchart(selected.event$KO, add=TRUE, at=1, pch=19, cex=1, vertical=TRUE, method="jitter", jitter=0.1)
stripchart(selected.event$WT, add=TRUE, at=2, pch=19, cex=1, vertical=TRUE, method="jitter", jitter=0.1)

selected.event <- list(KO=unlist(rmats.se.filtered["10192",c("PSI_S1R1","PSI_S1R2","PSI_S1R3")]), WT=unlist(rmats.se.filtered["10187",c("PSI_S2R1","PSI_S2R2","PSI_S2R3")]))
boxplot(selected.event, col=c("red","cyan"), ylab="deltaPSI", main="Lef1", las=1, ylim=c(0,1))
set.seed(1234)
stripchart(selected.event$KO, add=TRUE, at=1, pch=19, cex=1, vertical=TRUE, method="jitter", jitter=0.1)
stripchart(selected.event$WT, add=TRUE, at=2, pch=19, cex=1, vertical=TRUE, method="jitter", jitter=0.1)

