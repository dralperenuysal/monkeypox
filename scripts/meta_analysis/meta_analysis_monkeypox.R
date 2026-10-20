library(DExMA)
library(Biobase)
library(glue)
library(limma)
library(edgeR)
library(impute)
library(dplyr)
library(tibble)
library(annotables)

# Define the series
# GSE22121 (vaccination) and GSE24943 (Ebola) are out of scope; GSE11234 holds the
# same 124 hybridizations as GSE24125, so only GSE24125 is kept.
series <- c("GSE36854", "GSE24125")

# Initialize an empty list to store ExpressionSet objects
esets <- list()
pheno.list <- list()

# Loop through each series in the 'series' vector
for (serie in series) {
  expression_file <- glue("output/exp_fixed/{serie}_exp_fixed.tsv")
  pheno_file <- glue("output/exp_fixed/{serie}_phenotype_fixed.tsv")
  
  # Check if both files exist
  if (file.exists(expression_file) && file.exists(pheno_file)) {
    # Use tryCatch to handle potential errors
    tryCatch({
      # Read expression data
      exprs <- read.delim(expression_file, sep = '\t', check.names = F) %>% distinct(symbol, .keep_all = T) # Ensure row names are set
      phenodata <- read.delim(pheno_file, sep = '\t', row.names = 1)
      exprs <- exprs[, rownames(phenodata)]
      exprs <- as.matrix(exprs)  # Convert to matrix
      
      # Read phenotype data
      pheno.list[[serie]] <- phenodata
      phenodata.adf <- AnnotatedDataFrame(phenodata)  # Convert to AnnotatedDataFrame
      
      exprs.imp <- impute.knn(exprs, k = 10)
      exprs <- exprs.imp$data
      
      # Create ExpressionSet and add to list
      esets[[serie]] <- ExpressionSet(assayData = exprs,
                                      phenoData = phenodata.adf)
    }, error = function(e) {
      # Handle the error: print a message and continue
      message(glue("Error processing series {serie}: {e$message}"))
    })
  } else {
    warning(glue("Files for series {serie} do not exist: {expression_file}, {pheno_file}"))
  }
}

meta.obj <- list()

# Loop through each eset in esets using names
for (name in names(esets)) {
  eset <- esets[[name]]
  
  meta.element <- elementObjectMA(expressionMatrix = eset,
                                  groupPheno = "state",
                                  expGroup = "MP",
                                  refGroup = "control")
  #imputation <- missGenesImput(meta.element, k=7)
  #meta.element.imput <- imputation$meta.element
  meta.obj[[name]] <- meta.element
}

meta.obj <- allSameID(meta.obj, finalID="GeneSymbol", organism = "Homo sapiens")
meta.obj <- dataLog(meta.obj)

png("output/plots/heterogenity.png", width = 7, height = 7, units = "in", res = 300)
print(heterogeneityTest(meta.obj))
dev.off()
effect.size <- as.data.frame(DExMA::calculateES(meta.obj))

resultsMA <- metaAnalysisDE(meta.obj, typeMethod="REM")
results.df <- as.data.frame(resultsMA)
results.df.filtered <- results.df %>%
  filter(propDataset == 1) %>%
  filter(Zval > 2 | Zval < -2) %>%
  filter(FDR < 0.05)

# Filter the protein_coding genes
annot <- grch38 %>% dplyr::select(c(3,4,5,6,8,9)) %>%
  filter(biotype == "protein_coding")

results.df.filtered$symbol <- rownames(results.df.filtered)
results.df.filtered <- results.df.filtered %>%
  filter(symbol %in% annot$symbol)
results.merged <- merge(results.df.filtered, annot)
results.merged <- results.merged %>% distinct(symbol, .keep_all = T)


# Heatmap only from protein-coding DEGs
png("output/plots/heatmap_10_down.png", width = 9, height = 6, units = "in", res = 300)
makeHeatmap(objectMA =  meta.obj,
            numSig = 10,
            logFCSig = 0,
            fdrSig = 0.05,
            show_rownames = T,
            resMA = results.merged %>% column_to_rownames("symbol"),
            regulation = "down",
            scaling = "zscor")
dev.off()

downregulated <- filter(results.merged, AveFC < 0)$symbol
upregulated <- filter(results.merged, AveFC > 0)$symbol
# Enrichment background: protein-coding genes tested in all datasets
universe <- intersect(rownames(filter(results.df, propDataset == 1)), annot$symbol)

write.table(downregulated, file = "output/ppi_hub/downregulated.txt", row.names = F, col.names = F, quote = F)
write.table(upregulated, file = "output/ppi_hub/upregulated.txt", row.names = F, col.names = F, quote = F)
write.table(universe, file = "output/results/universe.txt", row.names = F, col.names = F, quote = F)
write.table(results.merged, file = "output/results/merged_results.tsv", sep = '\t')