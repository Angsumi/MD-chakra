# run_gsea_human_orthologs.R
# Script to perform GSEA on a non-model organism by mapping to Human orthologs

dir.create("/home/angsuman/extra_spac/R_libs", recursive=TRUE, showWarnings=FALSE)
.libPaths(c("/home/angsuman/extra_spac/R_libs", .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos="http://cran.us.r-project.org", lib="/home/angsuman/extra_spac/R_libs")
}
packages <- c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, lib="/home/angsuman/extra_spac/R_libs", ask=FALSE, update=FALSE)
  }
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

cat("Loading DESeq2 results...\n")
res <- read.csv("downstream_results/10_Annotated_DESeq2_results.csv", row.names=1)

# Extract Gene Symbols from the GTF annotation if they exist
# In our DESeq2 results, the row names are gene_ids (e.g. NPIL_101111)
# We can parse the GTF to create a mapping of gene_id -> gene_name
cat("Parsing GTF for gene symbols...\n")
gtf <- readLines("annotation.gtf")
gtf_genes <- gtf[grepl("\tgene\t", gtf)]

gene_map <- data.frame(gene_id=character(), gene_name=character(), stringsAsFactors=FALSE)
# Simple parsing (could be optimized, but works for standard NCBI GTFs)
for(line in gtf_genes) {
  if(grepl("gene_id", line) && grepl("gene ", line)) {
    g_id <- sub('.*gene_id "([^"]+)".*', '\\1', line)
    g_name <- sub('.*gene "([^"]+)".*', '\\1', line)
    gene_map <- rbind(gene_map, data.frame(gene_id=g_id, gene_name=g_name))
  }
}
# Remove duplicates
gene_map <- unique(gene_map)

# Merge the gene names into our results
res$gene_id <- rownames(res)
res <- merge(res, gene_map, by="gene_id", all.x=TRUE)

# Drop rows with no gene name
res <- res[!is.na(res$gene_name) & res$gene_name != "", ]

# Convert spider gene names to HUMAN gene symbols (UPPERCASE)
# This acts as our ortholog mapping!
res$human_symbol <- toupper(res$gene_name)

# We need a ranked list of genes for GSEA. We will rank by log2FoldChange.
# Remove NAs in log2FoldChange
res <- res[!is.na(res$log2FoldChange), ]
res <- res[order(res$log2FoldChange, decreasing = TRUE), ]

# Create a named vector for GSEA
geneList <- res$log2FoldChange
names(geneList) <- res$human_symbol

# Remove duplicate symbols (keep the one with the highest absolute fold change)
geneList <- geneList[!duplicated(names(geneList))]

cat("Starting GSEA using Human GO terms...\n")
gse <- gseGO(geneList=geneList, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

outdir <- "downstream_results/GSEA"
dir.create(outdir, showWarnings=FALSE)

if (nrow(as.data.frame(gse)) > 0) {
  cat("GSEA found significant enriched terms!\n")
  write.csv(as.data.frame(gse), file.path(outdir, "14_GSEA_Enrichment_Results.csv"))
  
  # Generate plots
  library(ggplot2)
  p1 <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  ggsave(file.path(outdir, "14_GSEA_Dotplot.pdf"), plot=p1, width=8, height=6)
  
  cat("GSEA plots saved to downstream_results/GSEA\n")
} else {
  cat("No significant enrichment found for GSEA at the current threshold.\n")
  cat("This is expected since the initial DESeq2 run found 0 differentially expressed genes.\n")
}
