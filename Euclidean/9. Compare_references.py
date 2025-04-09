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

# Load mobility
path = "/storage/enyaa/REVISED/mobility_classification_all_wrong.csv"
mobility_df = pd.read_csv(path, sep=",")

# Filter for the chromosomal genes
chromosomal_df = mobility_df[mobility_df['Gene_name'].isin(chromosomal_genes)]

# Keep only the "not mobile" genes
not_mobile_df = chromosomal_df[chromosomal_df['Mobility'] == 'Not_mobile']


# Remove SME-2 and QnrC (not chromosomal)
genes_to_remove = ["SME-2", "QnrC"]
not_mobile_df = not_mobile_df.loc[~not_mobile_df["Gene_name"].isin(genes_to_remove)]

# Loop through taxonomy results
taxonomy_df_list = []
for gene in not_mobile_df['Gene_name']:
    #gene_name = 
    file_path = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/taxonomy_results_{gene}.csv"
    if os.path.exists(file_path):
        taxonomy_df = pd.read_csv(file_path, sep=",", header=None) 
        taxonomy_df.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        taxonomy_df = taxonomy_df.drop_duplicates()
        taxonomy_df_list.append(taxonomy_df)
    else:
        continue    
    
taxonomy_df_all = pd.concat(taxonomy_df_list)

not_mobile_taxonomy = not_mobile_df[['Gene_name']].merge(taxonomy_df_all, on='Gene_name', how='left') #Take 
gene_counts = not_mobile_taxonomy["Gene_name"].value_counts()
gene_counts.to_csv("/home/jolunds/newtest/compatible_genes")

keep_gene_counts = gene_counts[gene_counts >= 5].index # Remove genes with count lower than 5
compatible_df = not_mobile_taxonomy[not_mobile_taxonomy["Gene_name"].isin(keep_gene_counts)]

# Load euclidean df
euclidean_df_list = []
for gene in compatible_df['Gene_name'].unique():
    file_path = f"/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_{gene}.pkl"
    gene_euclidean_df = pd.read_pickle(file_path) #.T.reset_index()
    gene_euclidean_df.insert(0, 'Gene_name', gene)# Add gene name column
    gene_euclidean_df.rename(columns={'index': 'Bacteria_ID', f'{gene}': 'Euclidean_distance'}, inplace=True)

    bacteria_ids = compatible_df[compatible_df['Gene_name'] == gene]['Bacteria_ID']
    gene_euclidean_df = gene_euclidean_df[gene_euclidean_df['Bacteria_ID'].isin(bacteria_ids)]
    euclidean_df_list.append(gene_euclidean_df)

comp_euclidean_df = pd.concat(euclidean_df_list).reset_index(drop=True)
comp_euclidean_df['Reference'] = 'Compatible' # Add 'Compatible'



# INCOMPATIBLE ------------------------
gram_positive_genera = ["Clostridium", "Bacillus", "Clostridioides", "Listeria", "Staphylococcus", "Enterococcus", "Lactobacillus", "Leuconostoc",
                        "Streptococcus"]

# Load full taxonomy and filter for gram_positive genera
path = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_lineage_df = pd.read_csv(path, sep="\t", header=None)
full_lineage_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

genus_mapping = full_lineage_df[["Bacteria_ID", "Genus"]]

gram_positive_bacteria = full_lineage_df[full_lineage_df["Genus"].isin(gram_positive_genera)]
incomp_bacteria_id = list(gram_positive_bacteria["Bacteria_ID"])

random.seed(42)
incomp_bacteria_id = random.sample(incomp_bacteria_id, k=10_000)
incomp_gene_names = ["NDM", "IMP", "GIM", "SPM", "VIM"] #Metallo betalaktamases 

# Ladda gene namn
with open("/storage/jolunds/REVISED/gene_names.txt", "r") as f:
    all_genes = [line.strip() for line in f]

incomp_genes = [gene for gene in all_genes if any(gene.startswith(prefix) for prefix in incomp_gene_names)]
random.seed(50)
incomp_genes = random.sample(incomp_genes, k=75)

# Loopa gennamn, filtera på bacteria id incomp_bacteria_id
incomp_euclidean_list = []
for gene in incomp_genes:
    file_path = f"/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_{gene}.pkl"
    gene_euclidean_df = pd.read_pickle(file_path) #.T.reset_index()
    gene_euclidean_df.insert(0, 'Gene_name', gene)# Add gene name column
    gene_euclidean_df.rename(columns={'index': 'Bacteria_ID', f'{gene}': 'Euclidean_distance'}, inplace=True)
    
    gene_euclidean_df = gene_euclidean_df[gene_euclidean_df['Bacteria_ID'].isin(incomp_bacteria_id)]
    
    #append
    incomp_euclidean_list.append(gene_euclidean_df)
    
incomp_euclidean_df = pd.concat(incomp_euclidean_list).reset_index(drop=True)
incomp_euclidean_df['Reference'] = 'Incompatible'

genes = incomp_euclidean_df['Gene_name'].nunique()
print(incomp_euclidean_df.head(10))
print(len(incomp_euclidean_df))
print(genes)

# Random 12k from incomp euclidean
sample_incomp_euclidean_df = incomp_euclidean_df.sample(n=500, random_state=42)

# Wilcoxon rank sum test
group1 = list(comp_euclidean_df['Euclidean_distance'])
group2 = list(sample_incomp_euclidean_df['Euclidean_distance'])

# Perform a one-sided Mann-Whitney U test (alternative='less' tests if group1 < group2)
stat, p_value = mannwhitneyu(group1, group2, alternative='less')  # Test if group1 has smaller values

print(f"Statistic: {stat}, P-value: {p_value}")

# Combine incomp and comp
reference_euclidean_df = pd.concat([comp_euclidean_df, sample_incomp_euclidean_df], ignore_index=True)
reference_euclidean_df = reference_euclidean_df.sort_values(by='Euclidean_distance').reset_index(drop=True)
reference_euclidean_df = reference_euclidean_df.merge(genus_mapping, on=['Bacteria_ID'], how='left')
reference_euclidean_df.to_csv("/home/jolunds/newtest/reference_euclidean_df.csv")
# Plot results
plt.figure(figsize=(10, 6))
sns.histplot(data=reference_euclidean_df, x='Euclidean_distance', hue='Reference', multiple='stack', bins=30)

plt.xlabel("Euclidean distance")
plt.ylabel("Number of bacteria")
plt.title(f"p={p_value}")

#plt.savefig("/home/enyaa/gene_genome/histogram_references.png")
plt.savefig("/home/jolunds/newtest/histogram_references_filtered.png")
plt.close()