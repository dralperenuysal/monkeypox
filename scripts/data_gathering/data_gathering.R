library(GEOquery)
library(DExMA)

# Read the series
series <- read.delim("resource/series.txt", header = F, col.names = "serie")
series <- series$serie

# Define the empty lists
sample.informations <- list()
expression.sets <- list()
exps <- list()
platforms <- list()

# Get the sample information 
for (serie in series) {
  gse <- getGEO(GEO = serie, GSEMatrix = T, AnnotGPL = T)
  gse <- gse[[1]]
  
  expression.sets[serie] <- gse
  exps[[serie]] <- exprs(gse)
  
  sample.info <- pData(phenoData(gse))
  sample.informations[[serie]] <- sample.info
  
  pl <- gse@annotation
  plat <- getGEO(pl)
  plat <- plat@dataTable@table
  platforms[[serie]] <- plat
}

# Save the sample informations
for (d in c("phenotype", "platforms", "expressions")) dir.create(file.path("resource", d), recursive = TRUE, showWarnings = FALSE)
for (x in names(sample.informations)) {
  write.table(sample.informations[[x]], file = paste0("resource/phenotype/", x, "_phenotype_raw.tsv"), sep = '\t')
  write.table(platforms[[x]], file = paste0("resource/platforms/", x, "_platform.tsv"), sep = '\t')
  write.table(exps[[x]], file = paste0("resource/expressions/", x, "_expression_raw.tsv"), sep = '\t')
}
