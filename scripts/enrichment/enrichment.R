library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(dplyr)

to_entrez <- function(file) {
  genes <- read.delim(file, header = F, col.names = "gene")$gene
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
}

# All downregulated DEGs; genome-wide background (clusterProfiler default)
genes <- to_entrez("output/ppi_hub/downregulated.txt")
universe <- to_entrez("output/results/universe.txt")
cat("Input genes:", length(genes), " Array universe:", length(universe), "\n")

save_enrichment <- function(res, name) {
  out <- res@result %>% filter(p.adjust < 0.05) %>%
    mutate(p.adjust = signif(p.adjust, 3), qvalue = signif(qvalue, 3))
  write.table(out, file = paste0("output/enrichment/", name, "_enrichment.tsv"), sep = '\t', quote = F, row.names = F)
  cat(name, ":", nrow(out), "significant terms\n")
  if (nrow(out) < 2) return(invisible())
  ggsave(paste0("output/enrichment/", name, "_dotplot.png"), dotplot(res, showCategory = 15), width = 10, height = 8, dpi = 300)
  ggsave(paste0("output/enrichment/", name, "_emapplot.png"), emapplot(pairwise_termsim(res), showCategory = 30),
         width = 10, height = 8, dpi = 300)
}

for (ont in c("BP", "CC", "MF")) {
  go <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = ont,
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  save_enrichment(go, ont)
}

kegg <- enrichKEGG(gene = genes, organism = 'hsa',
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
save_enrichment(setReadable(kegg, org.Hs.eg.db, keyType = "ENTREZID"), "kegg")

# Sensitivity: background restricted to genes measured in all datasets
for (ont in c("BP", "CC", "MF")) {
  go <- enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db, ont = ont, pAdjustMethod = "BH")
  cat("Array background", ont, ":", sum(go@result$p.adjust < 0.05), "terms with FDR < 0.05\n")
}
kegg <- enrichKEGG(gene = genes, universe = universe, organism = 'hsa', pAdjustMethod = "BH")
cat("Array background KEGG :", sum(kegg@result$p.adjust < 0.05), "terms with FDR < 0.05\n")
