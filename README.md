# Meta-analysis of host gene expression in monkeypox virus infection

Code and results for the manuscript *Integrated Meta-analysis and Functional Enrichment of Monkeypox
Virus-infection Gene Expression Datasets: Identification of Key Pathways and Hub Genes*.

Random-effects meta-analysis of two GEO microarray datasets of human cells infected with monkeypox virus (MPXV)
in vitro, followed by STRING protein–protein interaction network analysis and GO/KEGG enrichment.

## Datasets

| GEO series | Cells | MPXV-infected | Uninfected controls |
|---|---|---|---|
| [GSE24125](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24125) | Primary human macrophages and dermal fibroblasts, 0–48 h | 24 | 28 |
| [GSE36854](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36854) | HeLa | 2 | 2 |

Excluded series: GSE22121 (vaccination samples only), GSE24943 (no MPXV-infected samples) and GSE11234 (same
hybridizations as GSE24125). Within GSE24125, samples infected with other viruses or exposed to inactivated MPXV
are excluded.

## Pipeline

Run all steps from the repository root.

1. `Rscript scripts/data_gathering/data_gathering.R`: download series matrices and platform annotations
   (series listed in `resource/series.txt`) into `resource/`.
2. `python scripts/phenotype_fix/<GSE>_phenotype_fix.py` and `python scripts/expression_fix/<GSE>_expression_fix.py`
   for GSE24125 and GSE36854: select samples and map probes to gene symbols (`output/exp_fixed/`).
3. `Rscript scripts/meta_analysis/meta_analysis_monkeypox.R`: KNN imputation (k = 10), DExMA random-effects
   meta-analysis, heterogeneity Q-Q plot and heatmap (`output/results/`, `output/plots/`, DEG lists in `output/ppi_hub/`).
4. `Rscript scripts/ppi/ppi_hub.R`: STRING v12 network (score ≥ 0.900) and top 10 hub genes by betweenness
   centrality (`output/ppi_hub/`).
5. `Rscript scripts/enrichment/enrichment.R`: GO (BP, CC, MF) and KEGG enrichment of the downregulated genes,
   plus a sensitivity analysis with the array-restricted background (`output/enrichment/`).

Steps 4 and 5 download STRING and KEGG data and need internet access.

## Main results

- `output/results/merged_results.tsv`: 309 protein-coding DEGs (306 downregulated, 3 upregulated)
- `output/ppi_hub/hub_genes_10.csv`: hub genes
- `output/enrichment/*_enrichment.tsv`: significant GO terms and KEGG pathways

## Software

R 4.3.3 with GEOquery 2.70.0, Biobase 2.62.0, impute 1.76.0, DExMA 1.10.7, annotables 0.2.0, dplyr 1.1.4,
STRINGdb 2.14.3, igraph 2.1.3, clusterProfiler 4.10.1, org.Hs.eg.db 3.18.0, enrichplot 1.22.0 and DOSE 3.28.2;
Python with pandas 2.3.3 and NumPy 2.4.6.

## License

Code is released under the [MIT License](LICENSE). GEO data are subject to the terms of the original submitters.
