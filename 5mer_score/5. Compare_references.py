import json
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import euclidean
import random
import numpy as np
import os
from scipy.stats import mannwhitneyu

# COMPATIBLE --------------------------------------------
with open("/storage/jolunds/card.json", "r") as file:
    card_info = json.load(file)
search_word = "chromosom" # Search for chromosom

# Find names of chromosomal genes
chromosomal_genes = []
for key, value in card_info.items():
    if search_word in str(value):
        # Check if 'model_name' exists in the value 
        if isinstance(value, dict) and 'model_name' in value:
           chromosomal_genes.append(value['model_name'])


# Loop through chromosomal genes
not_mobile_genes = []
for gene_name in chromosomal_genes:
    file = f"{gene_name}_euclidean.tsv" #load euclidean distance file
    gene_euclidean_df = pd.read_csv(file, sep="\t")

    # Remove non-matches
    gene_euclidean_df = gene_euclidean_df[gene_euclidean_df["Match_status"] == "Match"]

    # Add gene name column
    gene_euclidean_df.insert(0, "Gene_name", gene_name)

    # Count how many unqiue for each taxonomy level
    taxonomy_levels = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    taxonomy_counts = gene_euclidean_df[taxonomy_levels].nunique() # group by Gene ID and count the unique values for different taxonomy levels, and count how many bacteria the gene is found in

    if taxonomy_counts["Genus"] == 1:
        not_mobile_genes.append(gene_euclidean_df) # Has to be in only one genus, but can be in multiple species  

# Concat
compatible_df = pd.concat(not_mobile_genes).reset_index(drop=True)
compatible_df['Reference'] = "Compatible"

# INCOMPATIBLE ------------------------

# Load gene names
with open("/storage/enyaa/FINAL/gene_names.txt", "r") as f:
    all_genes = [line.strip() for line in f]

# Define gene groups
gene_groups = {
    "Bacillota": ["NDM", "IMP", "VIM"],
    "Pseudomonadota": ["van"]
}

incompatible_list = []

for phylum, gene_prefixes in gene_groups.items():
    # Get genes that match any of the given prefixes
    gene_matches = [gene for gene in all_genes if any(gene.startswith(prefix) for prefix in gene_prefixes)]
    
    # Sample 75 genes for each group if needed (or adjust as desired)
    random.seed(50)
    gene_matches = random.sample(gene_matches, k=50 if len(gene_matches) >= 50 else len(gene_matches))
    
    for gene in gene_matches:
        # Load Euclidean data
        file_path = f"/storage/enyaa/FINAL/KMER/euclidean_split_genes_filtered/euclidean_df_{gene}.pkl"
        gene_euclidean_df = pd.read_pickle(file_path)
        gene_euclidean_df.insert(0, 'Gene_name', gene)
        
        # Append
        incompatible_list.append(gene_euclidean_df)

# Combine and label
incompatible_df = pd.concat(incompatible_list).reset_index(drop=True)
incompatible_df['Reference'] = 'Incompatible'

# Random 500 from incomp euclidean
shared_keys = incompatible_df[['Gene_name', 'Bacteria_ID']].drop_duplicates()
sample_keys = shared_keys.sample(n=500, random_state=42)

# Step 2: Filter both DataFrames using the sampled keys
sample_incompatible_df = incompatible_df.merge(sample_keys, on=['Gene_name', 'Bacteria_ID'], how='inner').drop_duplicates()

# Wilcoxon rank sum test
group1 = list(compatible_df['Euclidean_distance'])
group2 = list(sample_incompatible_df['Euclidean_distance'])

# Perform a one-sided Mann-Whitney U test (alternative='less' tests if group1 < group2)
stat, p_value = mannwhitneyu(group1, group2, alternative='less')  # Test if group1 has smaller values

# Combine incomp and comp
reference_df = pd.concat([compatible_df, sample_incompatible_df], ignore_index=True)
reference_df = reference_df.sort_values(by='Euclidean_distance').reset_index(drop=True)

# Plot results
plt.figure(figsize=(8, 6))
g = sns.histplot(data=reference_df, x='Euclidean_distance', hue='Reference', hue_order=["Compatible", "Incompatible"], multiple='stack', bins=30, palette={'Compatible': 'mediumseagreen', 'Incompatible': 'palevioletred'})

plt.xlabel("5mer score", fontsize=16)
plt.setp(g.get_legend().get_texts(), fontsize=16)
g.legend_.set_title(None)
plt.ylabel("", fontsize=16)
plt.tight_layout()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.savefig("/storage/jolunds/FINAL/REFERENCES/5mer_score_references.png")
plt.close()

print(f"p-value for 5mer score: {p_value}")
reference_df.to_csv("/storage/jolunds/FINAL/REFERENCES/5mer_score_reference_df.csv")
