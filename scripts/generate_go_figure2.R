# generate_go_figure2.R
dir.create("/home/angsuman/extra_spac/R_libs", recursive=TRUE, showWarnings=FALSE)
.libPaths(c("/home/angsuman/extra_spac/R_libs", .libPaths()))

suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

cat("Loading DESeq2 results...\n")
res <- read.csv("downstream_results/10_Annotated_DESeq2_results.csv", row.names=1)

cat("Parsing GTF for gene symbols...\n")
gtf <- readLines("annotation.gtf")
gtf_genes <- gtf[grepl("\tgene\t|product", gtf)]

gene_map <- data.frame(gene_id=character(), gene_name=character(), stringsAsFactors=FALSE)
for(line in gtf_genes) {
  if(grepl("gene_id", line) && (grepl("gene ", line) || grepl("product", line))) {
    g_id <- sub('.*gene_id "([^"]+)"*;.*', '\\1', line)
    g_name <- ""
    if(grepl("gene ", line)) {
      g_name <- sub('.*gene "([^"]+)".*', '\\1', line)
    } else {
      g_name <- sub('.*product "([^"]+)".*', '\\1', line)
    }
    # Keep first word as symbol
    g_name <- strsplit(g_name, " ")[[1]][1]
    gene_map <- rbind(gene_map, data.frame(gene_id=g_id, gene_name=g_name))
  }
}
gene_map <- unique(gene_map)

# Merge
res$gene_id <- rownames(res)
res <- merge(res, gene_map, by="gene_id", all.x=TRUE)
res <- res[!is.na(res$gene_name) & res$gene_name != "", ]
res$human_symbol <- toupper(res$gene_name)

# All expressed genes
all_genes <- unique(res$human_symbol)

# DEGs (Since padj < 0.05 gives 0, we use pvalue < 0.05 and abs(LFC) > 1 to show something)
degs <- res[!is.na(res$pvalue) & res$pvalue < 0.05 & abs(res$log2FoldChange) > 1, ]
deg_genes <- unique(degs$human_symbol)

cat("Converting to Entrez IDs...\n")
all_entrez <- bitr(all_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
deg_entrez <- bitr(deg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID

# Function to run groupGO
get_go_level2 <- function(genes, ont) {
  gg <- groupGO(gene=genes, OrgDb=org.Hs.eg.db, ont=ont, level=2, readable=TRUE)
  df <- as.data.frame(gg)
  df$Ontology <- ont
  return(df)
}

cat("Grouping GO Terms (Level 2)...\n")
all_bp <- get_go_level2(all_entrez, "BP")
all_cc <- get_go_level2(all_entrez, "CC")
all_mf <- get_go_level2(all_entrez, "MF")
all_go <- rbind(all_bp, all_cc, all_mf)
all_go$Type <- "All Unigene"

deg_bp <- get_go_level2(deg_entrez, "BP")
deg_cc <- get_go_level2(deg_entrez, "CC")
deg_mf <- get_go_level2(deg_entrez, "MF")
deg_go <- rbind(deg_bp, deg_cc, deg_mf)
deg_go$Type <- "DEG Unigene"

# Combine
final_df <- rbind(all_go, deg_go)
# Filter empty categories
final_df <- final_df[final_df$Count > 0, ]

# Select top categories to keep plot readable (top 15 per ontology based on All Unigene count)
top_terms <- final_df %>% 
  filter(Type == "All Unigene") %>% 
  group_by(Ontology) %>% 
  top_n(15, Count) %>% 
  pull(Description)

plot_df <- final_df[final_df$Description %in% top_terms, ]

# Plot
cat("Generating Figure 2 Plot...\n")
plot_df$Ontology <- factor(plot_df$Ontology, levels=c("BP", "CC", "MF"), labels=c("Biological Process", "Cellular Component", "Molecular Function"))

p <- ggplot(plot_df, aes(x=Description, y=Count, fill=Type)) +
  geom_bar(stat="identity", position="dodge") +
  facet_grid(~Ontology, scales="free_x", space="free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8),
        strip.background = element_rect(fill="gray90"),
        legend.position = "top",
        legend.title = element_blank()) +
  scale_fill_manual(values=c("All Unigene"="#4B8BBE", "DEG Unigene"="#E06666")) +
  labs(x="", y="Number of Genes", title="Gene Ontology Classification")

ggsave("Figure2_GO_Classification.pdf", plot=p, width=12, height=6)
cat("Saved to Figure2_GO_Classification.pdf\n")
