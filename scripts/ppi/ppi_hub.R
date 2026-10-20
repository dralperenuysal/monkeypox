library(STRINGdb)
library(igraph)

dir.create("resource/string", recursive = TRUE, showWarnings = FALSE)

# STRING v12, Homo sapiens, highest confidence (0.900)
string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 900,
                          network_type = "full", input_directory = "resource/string")

genes <- read.delim("output/ppi_hub/downregulated.txt", header = F, col.names = "gene")
mapped <- string_db$map(genes, "gene", removeUnmappedRows = TRUE)

edges <- string_db$get_interactions(mapped$STRING_id)
edges <- unique(edges[, c("from", "to", "combined_score")])
id2gene <- setNames(mapped$gene, mapped$STRING_id)
edges$from <- id2gene[edges$from]
edges$to <- id2gene[edges$to]

g <- simplify(graph_from_data_frame(edges, directed = FALSE), edge.attr.comb = "max")

# Hub genes: top 10 by betweenness centrality (same algorithm as cytoHubba)
btw <- sort(betweenness(g), decreasing = TRUE)
hubs <- data.frame(Rank = 1:10, Gene = names(btw)[1:10], Betweenness = round(btw[1:10], 2),
                   Degree = degree(g)[names(btw)[1:10]])

write.table(edges, "output/ppi_hub/ppi_network.csv", sep = ",", row.names = F, quote = F)
write.csv(hubs, "output/ppi_hub/hub_genes_10.csv", row.names = F)
write.table(hubs$Gene, "output/ppi_hub/downreg_hub_genes.txt", row.names = F, col.names = F, quote = F)
cat("Mapped:", nrow(mapped), "/", nrow(genes), " Nodes:", vcount(g), " Edges:", ecount(g), "\n")
print(hubs)

# Network figures (hub genes highlighted)
set.seed(42)
lay <- layout_with_fr(g)
is_hub <- V(g)$name %in% hubs$Gene
png("output/ppi_hub/ppi_network_09hc.png", width = 12, height = 12, units = "in", res = 300)
plot(g, layout = lay, vertex.size = ifelse(is_hub, 6, 3),
     vertex.color = ifelse(is_hub, "#D7301F", "#9ECAE1"), vertex.frame.color = NA,
     vertex.label = V(g)$name, vertex.label.cex = ifelse(is_hub, 0.9, 0.55),
     vertex.label.color = "black", vertex.label.family = "sans", edge.color = "grey75")
dev.off()

sub <- induced_subgraph(g, hubs$Gene)
png("output/ppi_hub/STRING network_Betweenness_top10.png", width = 7, height = 7, units = "in", res = 300)
plot(sub, layout = layout_in_circle(sub), vertex.size = 30, vertex.color = "#D7301F",
     vertex.frame.color = NA, vertex.label.color = "white", vertex.label.cex = 0.75,
     vertex.label.family = "sans", edge.color = "grey50", edge.width = 2)
dev.off()
