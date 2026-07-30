.libPaths(c("~/R/libs", .libPaths()))
suppressMessages(library(DESeq2))
counts_data <- read.table("counts/gene_counts.txt", header=TRUE, row.names=1, comment.char="#")
counts_data <- counts_data[, 6:9]
colnames(counts_data) <- c("NPFM1", "NPFM2", "NPFM3", "NPFM4")

coldata <- data.frame(row.names=colnames(counts_data), sample=colnames(counts_data))

cat("Creating DESeq2 object...\n")
dds <- DESeqDataSetFromMatrix(countData=counts_data, colData=coldata, design=~1)
dds <- estimateSizeFactors(dds)

cat("Running VST normalization...\n")
vsd <- vst(dds, blind=TRUE)
norm_counts <- assay(vsd)

write.csv(norm_counts, "vst_normalized_counts.csv")
cat("Saved to vst_normalized_counts.csv\n")
