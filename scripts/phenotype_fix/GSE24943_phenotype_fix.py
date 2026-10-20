import pandas as pd

"""NOT APPROPRIATE FOR THE ANALYSIS."""



# Define the GEO serie number
GEO = "GSE24943"

# Read the expression file
phenotype = pd.read_csv(f"resource/phenotype/{GEO}_phenotype_raw.tsv", sep="\t").reset_index()

# Comparsion column
state = "characteristics_ch2.1"

# Define the columns to be selected
cols = ['index', 'title', 'geo_accession', state]

# Choose the data
phenotype_chosen = phenotype[cols]

# Copy the targeted column
phenotype_chosen['state'] = phenotype_chosen[state]

# Filter the samples
phenotype_filtered = phenotype_chosen[phenotype_chosen['state'].str.contains("Monkeypox Virus|Mock infection")]

# Word to be searched
monkey = "Monkeypox"

# Filter and transform the dataset
phenotype_filtered['state'] = phenotype_filtered['state'].apply(lambda x: 'MP' if monkey in x else 'control')

# Rename the first column
phenotype_filtered = phenotype_filtered.rename(columns={'index': 'samples'})

# Save the expression data
phenotype_filtered.to_csv(f"output/exp_fixed/{GEO}_phenotype_fixed.tsv", sep='\t', index=False)