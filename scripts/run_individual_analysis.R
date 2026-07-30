# run_individual_analysis.R
# Comprehensive RNA-Seq downstream analysis script without pseudo-groups

# 1. Setup and load packages
.libPaths(c("~/R/libs", .libPaths()))

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

outdir <- "downstream_results"
dir.create(outdir, showWarnings=FALSE)

cat("Loading data...\n")
counts_file <- read.table("counts/gene_counts.txt", header=TRUE, row.names=1, comment.char="#")

# Extract annotation
gene_anno <- counts_file[, 1:5]
countData <- counts_file[, -(1:5)]
colnames(countData) <- gsub("alignments\\.|_sorted\\.bam", "", colnames(countData))

# Define experimental design with 4 independent samples (design = ~1)
colData <- data.frame(row.names = colnames(countData), sample = colnames(countData))

cat("Creating DESeq2 Dataset...\n")
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ 1)
mcols(dds)$basepairs <- gene_anno$Length

# Estimate size factors
dds <- estimateSizeFactors(dds)

cat("Calculating FPKM...\n")
fpkm_vals <- fpkm(dds)
write.csv(fpkm_vals, file.path(outdir, "6_FPKM_normalized_counts_individual.csv"))

# Transformation for PCA and Clustering
rld <- rlog(dds, blind=TRUE) # blind=TRUE is appropriate when there is no experimental design
vst_data <- assay(rld)
write.csv(vst_data, file.path(outdir, "normalized_counts_rlog.csv"))

cat("Generating PCA Plot...\n")
pcaData <- plotPCA(rld, intgroup=c("sample"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color=sample, label=sample)) +
  geom_point(size=4) +
  geom_text_repel() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA Plot (Individual Samples)") +
  theme_bw() +
  scale_color_brewer(palette="Set1")

ggsave(file.path(outdir, "7_PCA_plot_individual.pdf"), plot=p, width=6, height=5)

cat("Generating Sample Distance Plot...\n")
sampleDists <- dist(t(vst_data))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(sampleDistMatrix) <- colnames(vst_data)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf(file.path(outdir, "8_Sample_Distance_Matrix_individual.pdf"))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, main="Sample-to-Sample Distances")
dev.off()

cat("Generating Top 50 Variable Genes Heatmap...\n")
topVarGenes <- head(order(rowVars(vst_data), decreasing=TRUE), 50)
mat <- vst_data[ topVarGenes, ]
mat <- mat - rowMeans(mat)

# No annotation_col needed as samples are independent
pdf(file.path(outdir, "11_Heatmap_Top50_Variable_individual.pdf"))
pheatmap(mat, main="Top 50 Variable Genes", scale="row", show_colnames=TRUE)
dev.off()

cat("R analysis complete.\n")
