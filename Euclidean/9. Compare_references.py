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

'''
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

# Load mobility
path = "/storage/enyaa/REVISED/mobility_classification_all_wrong.csv"
mobility_df = pd.read_csv(path)

# Filter for the chromosomal genes
chromosomal_df = mobility_df[mobility_df['Gene_name'].isin(chromosomal_genes)]

# Keep only the "not mobile" genes
not_mobile_df = chromosomal_df[chromosomal_df['Mobility'] == 'Not_mobile']

# Remove SME-2 and QnrC (not chromosomal)
genes_to_remove = ["SME-2", "QnrC"]
not_mobile_df = not_mobile_df.loc[~not_mobile_df["Gene_name"].isin(genes_to_remove)]

not_mobile_df.to_csv("/storage/jolunds/compatible_genes.csv")
'''

# COMPATIBLE -------------------------
# Load compatible genes
compatible_genes_df = pd.read_csv("/storage/jolunds/compatible_genes.csv")

# Loop through taxonomy results
taxonomy_df_list = []
for gene in compatible_genes_df['Gene_name']:
    file_path = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/taxonomy_results_{gene}.csv"
    if os.path.exists(file_path):
        taxonomy_df = pd.read_csv(file_path, sep=",", header=None) 
        taxonomy_df.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        taxonomy_df = taxonomy_df.drop_duplicates()
        taxonomy_df_list.append(taxonomy_df)
    else:
        continue    
    
taxonomy_df_all = pd.concat(taxonomy_df_list)

compatible_df = compatible_genes_df[['Gene_name']].merge(taxonomy_df_all, on='Gene_name', how='left') #Take 
#gene_counts = not_mobile_taxonomy["Gene_name"].value_counts()

#keep_gene_counts = gene_counts[gene_counts >= 5].index # Remove genes with count lower than 5
#compatible_df = not_mobile_taxonomy[not_mobile_taxonomy["Gene_name"].isin(keep_gene_counts)]


# Load euclidean df and tRNA score
euclidean_df_list = []
tRNA_score_list = []
euclidean_500_list = []
for gene in compatible_df['Gene_name'].unique():
    file_path = f"/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_{gene}.pkl"
    gene_euclidean_df = pd.read_pickle(file_path) 
    gene_euclidean_df.insert(0, 'Gene_name', gene)# Add gene name column
    #gene_euclidean_df.rename(columns={'index': 'Bacteria_ID', f'{gene}': 'Euclidean_distance'}, inplace=True)

    bacteria_ids = compatible_df[compatible_df['Gene_name'] == gene]['Bacteria_ID']
    gene_euclidean_df = gene_euclidean_df[gene_euclidean_df['Bacteria_ID'].isin(bacteria_ids)]
    euclidean_df_list.append(gene_euclidean_df)
    
    tRNA_filepath = f"/storage/jolunds/REVISED/tRNA/tRNA_score_new/tRNA_score_{gene}.csv"
    gene_tRNA_df = pd.read_csv(tRNA_filepath)
    gene_tRNA_df.insert(0, 'Gene_name', gene) # Add gene name column
    gene_tRNA_df = gene_tRNA_df[gene_tRNA_df['Bacteria_ID'].isin(bacteria_ids)]
    tRNA_score_list.append(gene_tRNA_df)
    
    file_path_500 = f"/storage/enyaa/REVISED/KMER/FOR_GENE_LENGTH/euclidean_split_genes_500bp/euclidean_df_{gene}.pkl"
    if not os.path.exists(file_path_500):
        continue 
    
    gene_euclidean_500 = pd.read_pickle(file_path_500).T.reset_index()
    gene_euclidean_500 = gene_euclidean_500.rename(columns={'index': 'Bacteria_ID', gene: 'Euclidean_distance'})
    gene_euclidean_500.insert(0, 'Gene_name', gene)
    gene_euclidean_500 = gene_euclidean_500[gene_euclidean_500['Bacteria_ID'].isin(bacteria_ids)]
    euclidean_500_list.append(gene_euclidean_500)
    

comp_euclidean_df = pd.concat(euclidean_df_list).reset_index(drop=True)
comp_euclidean_df['Reference'] = 'Compatible' # Add 'Compatible'

comp_tRNA_df = pd.concat(tRNA_score_list).reset_index(drop=True)
comp_tRNA_df['Reference'] = 'Compatible' # Add 'Compatible'

comp_euclidean_500 = pd.concat(euclidean_500_list).reset_index(drop=True)
comp_euclidean_500['Reference'] = 'Compatible' # Add 'Compatible'


# INCOMPATIBLE ------------------------

# Load gene names
with open("/storage/jolunds/REVISED/gene_names.txt", "r") as f:
    all_genes = [line.strip() for line in f]

# Define gene groups
gene_groups = {
    "Bacillota": ["NDM", "IMP", "VIM"],
    "Pseudomonadota": ["van"]
}

incomp_euclidean_list = []
incomp_tRNA_list = []
incomp_euclidean_500_list = []
for phylum, gene_prefixes in gene_groups.items():
    # Get genes that match any of the given prefixes
    gene_matches = [gene for gene in all_genes if any(gene.startswith(prefix) for prefix in gene_prefixes)]
    
    # Sample 75 genes for each group if needed (or adjust as desired)
    random.seed(50)
    gene_matches = random.sample(gene_matches, k=50 if len(gene_matches) >= 50 else len(gene_matches))
    
    for gene in gene_matches:
        # Load Euclidean data
        file_path = f"/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_{gene}.pkl"
        gene_euclidean_df = pd.read_pickle(file_path)
        gene_euclidean_df.insert(0, 'Gene_name', gene)
        #gene_euclidean_df.rename(columns={'index': 'Bacteria_ID', f'{gene}': 'Euclidean_distance'}, inplace=True)
        
        # Load tRNA data
        tRNA_filepath = f"/storage/jolunds/REVISED/tRNA/tRNA_score_new/tRNA_score_{gene}.csv"
        gene_tRNA_df = pd.read_csv(tRNA_filepath)
        gene_tRNA_df.insert(0, 'Gene_name', gene)
        
        # Filter tRNA data to correct phylum
        gene_tRNA_df = gene_tRNA_df[gene_tRNA_df['Phylum'] == phylum]
        # Filter Euclidean data to matching bacteria IDs
        matching_ids = gene_tRNA_df['Bacteria_ID'].unique()
        gene_euclidean_df = gene_euclidean_df[gene_euclidean_df['Bacteria_ID'].isin(matching_ids)]
        
        # Append
        incomp_euclidean_list.append(gene_euclidean_df)
        incomp_tRNA_list.append(gene_tRNA_df)
        
        # Load euclidean 500 bp
        file_path_500 = f"/storage/enyaa/REVISED/KMER/FOR_GENE_LENGTH/euclidean_split_genes_500bp/euclidean_df_{gene}.pkl"
        if not os.path.exists(file_path_500):
            continue 
        gene_euclidean_500 = pd.read_pickle(file_path_500).T.reset_index()
        gene_euclidean_500 = gene_euclidean_500.rename(columns={'index': 'Bacteria_ID', gene: 'Euclidean_distance'})
        gene_euclidean_500.insert(0, 'Gene_name', gene)
        
        gene_euclidean_500 = gene_euclidean_500[gene_euclidean_500['Bacteria_ID'].isin(matching_ids)]
        incomp_euclidean_500_list.append(gene_euclidean_500)

# Combine and label
incomp_euclidean_df = pd.concat(incomp_euclidean_list).reset_index(drop=True)
incomp_euclidean_df['Reference'] = 'Incompatible'

incomp_tRNA_df = pd.concat(incomp_tRNA_list).reset_index(drop=True)
incomp_tRNA_df['Reference'] = 'Incompatible'

incomp_euclidean_500 = pd.concat(incomp_euclidean_500_list).reset_index(drop=True)
incomp_euclidean_500['Reference'] = 'Incompatible'

# Add Phylum info to Euclidean dataframe
incomp_euclidean_df = incomp_euclidean_df.merge(
    incomp_tRNA_df[['Bacteria_ID', 'Phylum']], 
    on='Bacteria_ID', 
    how='left'
)

incomp_euclidean_500 = incomp_euclidean_500.merge(
    incomp_tRNA_df[['Bacteria_ID', 'Phylum']], 
    on='Bacteria_ID', 
    how='left'
)

# Random 500 from incomp euclidean
shared_keys = incomp_euclidean_df[['Gene_name', 'Bacteria_ID']].drop_duplicates()
sample_keys = shared_keys.sample(n=500, random_state=42)

# Step 2: Filter both DataFrames using the sampled keys
sample_incomp_euclidean_df = incomp_euclidean_df.merge(sample_keys, on=['Gene_name', 'Bacteria_ID'], how='inner').drop_duplicates()
sample_incomp_tRNA_df = incomp_tRNA_df.merge(sample_keys, on=['Gene_name', 'Bacteria_ID'], how='inner')
sample_incomp_euclidean_500 = incomp_euclidean_500.merge(sample_keys, on=['Gene_name', 'Bacteria_ID'], how='inner').drop_duplicates()

# Wilcoxon rank sum test
group1_eu = list(comp_euclidean_df['Euclidean_distance'])
group2_eu = list(sample_incomp_euclidean_df['Euclidean_distance'])

tRNA_score = "tRNA_score_one_sided"
group1_tRNA = list(comp_tRNA_df[tRNA_score])
group2_tRNA = list(sample_incomp_tRNA_df[tRNA_score])

group1_500 = list(comp_euclidean_500['Euclidean_distance'])
group2_500 = list(sample_incomp_euclidean_500['Euclidean_distance'])

# Perform a one-sided Mann-Whitney U test (alternative='less' tests if group1 < group2)
stat_eu, p_value_eu = mannwhitneyu(group1_eu, group2_eu, alternative='less')  # Test if group1 has smaller values
stat_tRNA, p_value_tRNA = mannwhitneyu(group1_tRNA, group2_tRNA, alternative='less')  # Test if group1 has smaller values
stat_500, p_value_500 = mannwhitneyu(group1_500, group2_500, alternative='less')  # Test if group1 has smaller values

# Combine incomp and comp
reference_euclidean_df = pd.concat([comp_euclidean_df, sample_incomp_euclidean_df], ignore_index=True)
reference_euclidean_df = reference_euclidean_df.sort_values(by='Euclidean_distance').reset_index(drop=True)
reference_euclidean_df.to_csv("/home/jolunds/newtest/reference_euclidean_df.csv")

reference_tRNA_df = pd.concat([comp_tRNA_df, sample_incomp_tRNA_df], ignore_index=True)
reference_tRNA_df = reference_tRNA_df.sort_values(by=tRNA_score).reset_index(drop=True)
reference_tRNA_df.to_csv("/home/jolunds/newtest/reference_tRNA_df.csv")

reference_euclidean_500 = pd.concat([comp_euclidean_500, sample_incomp_euclidean_500], ignore_index=True)
reference_euclidean_500 = reference_euclidean_500.sort_values(by='Euclidean_distance').reset_index(drop=True)
reference_euclidean_500.to_csv("/home/jolunds/newtest/reference_euclidean_500.csv")

# Plot results
# Euclidean 
plt.figure(figsize=(10, 6))
g = sns.histplot(data=reference_euclidean_df, x='Euclidean_distance', hue='Reference', hue_order=["Compatible", "Incompatible"], multiple='stack', bins=30, palette={'Compatible': 'mediumseagreen', 'Incompatible': 'palevioletred'})

plt.xlabel("Euclidean distance")
plt.ylabel("Number of bacteria")
plt.title(f"p={p_value_eu}")

plt.savefig("/home/jolunds/newtest/References/euclidean_references.png")
plt.close()

# tRNA
plt.figure(figsize=(10, 6))
sns.histplot(data=reference_tRNA_df, x=tRNA_score, hue='Reference', hue_order=["Compatible", "Incompatible"], multiple='stack', bins=30, palette={'Compatible': 'mediumseagreen', 'Incompatible': 'palevioletred'})

plt.xlabel("tRNA score")
plt.ylabel("Number of bacteria")
plt.title(f"p={p_value_tRNA}")

plt.savefig("/home/jolunds/newtest/References/tRNA_references.png")
plt.close()

# Euclidean 500 bp
plt.figure(figsize=(10, 6))
g = sns.histplot(data=reference_euclidean_500, x='Euclidean_distance', hue='Reference', hue_order=["Compatible", "Incompatible"], multiple='stack', bins=30, palette={'Compatible': 'mediumseagreen', 'Incompatible': 'palevioletred'})

plt.xlabel("Euclidean distance 500 bp")
plt.ylabel("Number of bacteria")
plt.title(f"p={p_value_500}")

plt.savefig("/home/jolunds/newtest/References/euclidean_references_500.png")
plt.close()