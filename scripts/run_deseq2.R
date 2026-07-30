# run_deseq2.R
# Script for Differential Gene Expression Analysis using DESeq2

# Load required library
# Create and set user library for R packages
dir.create("~/R/libs", recursive=TRUE, showWarnings=FALSE)
.libPaths(c("~/R/libs", .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org", lib="~/R/libs")
if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2", lib="~/R/libs", ask=FALSE, update=FALSE)

library(DESeq2)

# 1. Load the featureCounts output
# featureCounts output has some metadata in the first 6 columns. 
# The count data starts from column 7.
countData <- read.table("counts/gene_counts.txt", header=TRUE, row.names=1, comment.char="#")

# Remove the metadata columns (Chr, Start, End, Strand, Length)
# Keep only the sample columns
countData <- countData[ , -c(1:5)]

# Clean up column names to match your sample names
# (e.g. changing 'alignments.NPFM1_sorted.bam' to 'NPFM1')
colnames(countData) <- gsub("alignments\\.|_sorted\\.bam", "", colnames(countData))
print("Sample columns found:")
print(colnames(countData))

# 2. Define your experimental design (Conditions)
# IMPORTANT: Update these based on your actual experiment!
# Example: Assuming NPFM1, NPFM2 are "Control", and NPFM3, NPFM4 are "Treatment"
sampleCondition <- factor(c("Control", "Control", "Treatment", "Treatment"))

# Create the colData dataframe required by DESeq2
colData <- data.frame(row.names = colnames(countData), condition = sampleCondition)
print("Experimental Design:")
print(colData)

# 3. Create DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

# Set the reference level for the condition (What are you comparing against?)
dds$condition <- relevel(dds$condition, ref = "Control")

# 4. Run the DESeq2 analysis pipeline
cat("Running DESeq2 pipeline...\n")
dds <- DESeq(dds)

# 5. Extract results
res <- results(dds)

# Order results by adjusted p-value (FDR)
resOrdered <- res[order(res$padj), ]

# Summary of the results
cat("\nSummary of Differential Expression (FDR < 0.1):\n")
summary(res)

# 6. Save results to a CSV file
write.csv(as.data.frame(resOrdered), file = "DESeq2_results.csv")
cat("\nResults saved to 'DESeq2_results.csv'\n")

# 7. (Optional) Basic Visualization - MA Plot
pdf("MA_plot.pdf")
plotMA(res, main="DESeq2: Treatment vs Control", ylim=c(-5,5))
dev.off()
cat("MA Plot saved to 'MA_plot.pdf'\n")
