# comprehensive_deseq2.R
# Comprehensive RNA-Seq downstream analysis script

# 1. Setup and load packages
dir.create("~/R/libs", recursive=TRUE, showWarnings=FALSE)
.libPaths(c("~/R/libs", .libPaths()))

packages <- c("DESeq2", "pheatmap", "RColorBrewer", "ggplot2", "ggrepel", "GenomicFeatures")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("pheatmap", "RColorBrewer", "ggplot2", "ggrepel")) {
      install.packages(pkg, repos="http://cran.us.r-project.org", lib="~/R/libs")
    } else {
      BiocManager::install(pkg, lib="~/R/libs", ask=FALSE, update=FALSE)
    }
  }
}

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

# Create output directory
outdir <- "downstream_results"
dir.create(outdir, showWarnings=FALSE)

# 2. Load featureCounts data and Annotation
cat("Loading data...\n")
counts_file <- read.table("counts/gene_counts.txt", header=TRUE, row.names=1, comment.char="#")

# Extract annotation (Chr, Start, End, Strand, Length) from featureCounts
gene_anno <- counts_file[, 1:5]
countData <- counts_file[, -(1:5)]
colnames(countData) <- gsub("alignments\\.|_sorted\\.bam", "", colnames(countData))

# Define experimental design
sampleCondition <- factor(c("Control", "Control", "Treatment", "Treatment"))
colData <- data.frame(row.names = colnames(countData), condition = sampleCondition)

# 3. Create DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "Control")

# To calculate FPKM we need the gene lengths. featureCounts provides it in column 'Length'
mcols(dds)$basepairs <- gene_anno$Length

# Run DESeq2
cat("Running DESeq2...\n")
dds <- DESeq(dds)
res <- results(dds)

# 4. Generate FPKM values (Requirement 6)
cat("Calculating FPKM...\n")
fpkm_vals <- fpkm(dds)
write.csv(fpkm_vals, file.path(outdir, "6_FPKM_normalized_counts.csv"))

# 5. Data Transformation for PCA and Clustering
rld <- rlog(dds, blind=FALSE)

# 6. PCA Plot (Requirement 7)
cat("Generating PCA Plot...\n")
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA Plot") +
  theme_bw()
ggsave(file.path(outdir, "7_PCA_plot.pdf"), plot=p, width=6, height=5)

# 7. Sample to Sample distance cluster plot (Requirement 8)
cat("Generating Sample Distance Plot...\n")
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rownames(sampleDistMatrix), sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(file.path(outdir, "8_Sample_Distance_Matrix.pdf"))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, main="Sample-to-Sample Distances")
dev.off()

# 8. Annotated DESeq2 Results (Requirement 10)
cat("Annotating Results...\n")
# Combine counts, annotations, and DESeq2 stats
res_df <- as.data.frame(res)
# Create linear fold change
res_df$FoldChange <- ifelse(res_df$log2FoldChange > 0, 2^res_df$log2FoldChange, -1/(2^res_df$log2FoldChange))

annotated_res <- cbind(gene_anno, countData, res_df)
write.csv(annotated_res, file.path(outdir, "10_Annotated_DESeq2_results.csv"))

# 9. Filter Cutoffs (Requirement 10 - 6 cutoffs)
cat("Filtering results by cutoffs...\n")
filters <- list(
  cutoff_1 = subset(annotated_res, padj < 0.05 & abs(log2FoldChange) >= 2),
  cutoff_2 = subset(annotated_res, padj < 0.05 & abs(log2FoldChange) >= 1.5),
  cutoff_3 = subset(annotated_res, padj < 0.05 & abs(log2FoldChange) >= 1),
  cutoff_4 = subset(annotated_res, pvalue < 0.05 & abs(log2FoldChange) >= 2),
  cutoff_5 = subset(annotated_res, pvalue < 0.05 & abs(log2FoldChange) >= 1.5),
  cutoff_6 = subset(annotated_res, pvalue < 0.05 & abs(log2FoldChange) >= 1)
)

for (name in names(filters)) {
  write.csv(filters[[name]], file.path(outdir, paste0("10_", name, "_filtered.csv")))
}

# 10. Complete set of up/down regulated genes (Requirement 9)
# Assuming cutoff 3 (FDR < 0.05, |LFC| >= 1) as standard for the complete set if not specified
upregulated <- subset(annotated_res, padj < 0.05 & log2FoldChange >= 1)
downregulated <- subset(annotated_res, padj < 0.05 & log2FoldChange <= -1)
write.csv(upregulated, file.path(outdir, "9_Upregulated_genes.csv"))
write.csv(downregulated, file.path(outdir, "9_Downregulated_genes.csv"))

# 11. Heatmap of 50 most variable genes (Requirement 11)
cat("Generating Top 50 Variable Genes Heatmap...\n")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), 50)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,"condition", drop=FALSE])
pdf(file.path(outdir, "11_Heatmap_Top50_Variable.pdf"))
pheatmap(mat, annotation_col=df, main="Top 50 Variable Genes", scale="row")
dev.off()

# 12. MA Plot (Requirement 11 - assuming 12 due to numbering)
cat("Generating MA Plot...\n")
pdf(file.path(outdir, "12_MA_Plot.pdf"))
plotMA(res, main="MA Plot (Treatment vs Control)", ylim=c(-5,5))
dev.off()

# 13. Volcano Plot (Requirement 12)
cat("Generating Volcano Plot...\n")
volcano_data <- as.data.frame(res)
volcano_data$gene <- rownames(volcano_data)
# Add significance labels
volcano_data$Significant <- "Not Significant"
volcano_data$Significant[which(volcano_data$padj < 0.05 & volcano_data$log2FoldChange >= 1)] <- "Upregulated"
volcano_data$Significant[which(volcano_data$padj < 0.05 & volcano_data$log2FoldChange <= -1)] <- "Downregulated"

vplot <- ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(pvalue), col=Significant)) +
  geom_point(alpha=0.4, size=1.5) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  theme_minimal() +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="black") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  ggtitle("Volcano Plot")
ggsave(file.path(outdir, "13_Volcano_Plot.pdf"), plot=vplot, width=6, height=5)

cat("Comprehensive analysis complete. Check the 'downstream_results' directory!\n")
