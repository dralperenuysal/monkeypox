import pandas as pd

# Define the GEO serie number
GEO = "GSE22121"

# Read the expression file
exp = pd.read_csv(f"resource/expressions/{GEO}_expression_raw.tsv", sep="\t").reset_index()

# Read the platform file
platform = pd.read_csv(f"resource/platforms/{GEO}_platform.tsv", sep="\t")

# Define the symbol column
symbol_column = "GB_ACC"

# Rename the first column
exp = exp.rename(columns={"index": "ID"})

platform_filtered = platform.dropna(subset=[symbol_column])
platform_filtered = platform_filtered[["ID", symbol_column]]

# Merge the platform information and the expression data
merged = pd.merge(exp, platform_filtered, on="ID")
merged = merged[merged.columns[1:]]
merged = merged.rename(columns={merged.columns[-1]: "symbol"})
merged = merged.fillna(0)

# Save the expression data
merged.to_csv(f"output/exp_fixed/{GEO}_exp_fixed.tsv", sep='\t', index=False)