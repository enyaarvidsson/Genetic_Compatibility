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
print(len(keep_gene_counts))
compatible_df = not_mobile_taxonomy[not_mobile_taxonomy["Gene_name"].isin(keep_gene_counts)]
#print(len(compatible_df))


# Load euclidean df
euclidean_df_list = []
for gene in compatible_df['Gene_name'].unique():
    file_path = f"/storage/enyaa/REVISED/KMER/euclidean_split_genes/euclidean_df_{gene}.pkl"
    gene_euclidean_df = pd.read_pickle(file_path).T.reset_index()
    gene_euclidean_df.insert(0, 'Gene_name', gene)# Add gene name column
    gene_euclidean_df.rename(columns={'index': 'Bacteria_ID', f'{gene}': 'Euclidean_distance'}, inplace=True)

    bacteria_ids = compatible_df[compatible_df['Gene_name'] == gene]['Bacteria_ID']
    gene_euclidean_df = gene_euclidean_df[gene_euclidean_df['Bacteria_ID'].isin(bacteria_ids)]
    euclidean_df_list.append(gene_euclidean_df)

comp_euclidean_df = pd.concat(euclidean_df_list).reset_index(drop=True)
comp_euclidean_df['Reference'] = 'Compatible' # Add 'Compatible'

print(comp_euclidean_df.head())
print(len(comp_euclidean_df))

# INCOMPATIBLE ------------------------
gram_positive_genera = ["Clostridium", "Bacillus", "Clostridioides", "Listeria", "Staphylococcus", "Enterococcus", "Lactobacillus", "Leuconostoc",
                        "Streptococcus"]

# Load full taxonomy and filter for gram_positive genera
path = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_lineage_df = pd.read_csv(path, sep="\t", header=None)
full_lineage_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

gram_positive_bacteria = full_lineage_df[full_lineage_df["Genus"].isin(gram_positive_genera)]
incomp_bacteria_id = list(gram_positive_bacteria["Bacteria_ID"])

random.seed(42)
incomp_bacteria_id = random.sample(incomp_bacteria_id, k=700)
incomp_gene_names = ["NDM", "IMP", "GIM", "SPM", "VIM"] #Metallo betalaktamases 

# Ladda gene namn
with open("/storage/jolunds/REVISED/gene_names.txt", "r") as f:
    all_genes = [line.strip() for line in f]

incomp_genes = [gene for gene in all_genes if any(gene.startswith(prefix) for prefix in incomp_gene_names)]
random.seed(50)
incomp_genes = random.sample(incomp_genes, k=120)
print(len(incomp_genes))

# Loopa gennamn, filtera på bacteria id incomp_bacteria_id
incomp_euclidean_list = []
for gene in incomp_genes:
    file_path = f"/storage/enyaa/REVISED/KMER/euclidean_split_genes/euclidean_df_{gene}.pkl"
    gene_euclidean_df = pd.read_pickle(file_path).T.reset_index()
    gene_euclidean_df.insert(0, 'Gene_name', gene)# Add gene name column
    gene_euclidean_df.rename(columns={'index': 'Bacteria_ID', f'{gene}': 'Euclidean_distance'}, inplace=True)
    
    gene_euclidean_df = gene_euclidean_df[gene_euclidean_df['Bacteria_ID'].isin(incomp_bacteria_id)]
    
    #append
    incomp_euclidean_list.append(gene_euclidean_df)
    
incomp_euclidean_df = pd.concat(incomp_euclidean_list).reset_index(drop=True)
incomp_euclidean_df['Reference'] = 'Incompatible'

print(incomp_euclidean_df.head())
print(len(incomp_euclidean_df))

###### WILCOXON

group1 = list(comp_euclidean_df['Euclidean_distance'])
group2 = list(incomp_euclidean_df['Euclidean_distance'])

print(group1[:10])
print(group2[:10])

# Perform a one-sided Mann-Whitney U test (alternative='less' tests if group1 < group2)
stat, p_value = mannwhitneyu(group1, group2, alternative='less')  # Test if group1 has smaller values

print(f"Statistic: {stat}, P-value: {p_value}")

# Combine incomp and comp
reference_euclidean_df = pd.concat([comp_euclidean_df, incomp_euclidean_df], ignore_index=True)
reference_euclidean_df = reference_euclidean_df.sort_values(by='Euclidean_distance').reset_index(drop=True)

reference_euclidean_df.to_csv("/home/jolunds/newtest/reference_euclidean_df.csv")
# Plot results
plt.figure(figsize=(10, 6))
sns.histplot(data=reference_euclidean_df, x='Euclidean_distance', hue='Reference', multiple='stack', bins=30)

plt.xlabel("Euclidean distance")
plt.ylabel("Number of bacteria")
plt.title(f"p={p_value}")

#plt.savefig("/home/enyaa/gene_genome/histogram_references.png")
plt.savefig("/home/jolunds/newtest/histogram_references.png")
plt.close()

'''
# Ensure only existing bacteria IDs are used
existing_bacteria_ids = set(euclidean_df.columns)  # Get available bacteria IDs
filtered_incomp_bacteria_id = incomp_bacteria_id[incomp_bacteria_id.isin(existing_bacteria_ids)]  # Filter valid ones

incomp_euclidean_df = euclidean_df.loc[:, euclidean_df.columns.isin(filtered_incomp_bacteria_id)]
incomp_euclidean_df = incomp_euclidean_df[incomp_euclidean_df.index.str.startswith(tuple(incomp_gene_names))]

# Merge to one 
comp_values = comp_euclidean_df.values.flatten()
incomp_values = incomp_euclidean_df.values.flatten()

df_plot = pd.DataFrame({
    'Euclidean Distance': list(comp_values) + list(incomp_values),
    'Reference': ['Compatible'] * len(comp_values) + ['Incompatible'] * len(incomp_values)
})


'''
'''
# Load gene dictionary and change to dataframe
with open("/storage/enyaa/REVISED/KMER/gene_dist/gene_kmer_distributions.pkl", "rb") as file: # "rb": read binary
    gene_dictionary = pickle.load(file) # for each gene name, we have the kmers (ex AAAAC) and the corresponding values (ex 0.0043) 

genes_df = pd.DataFrame.from_dict(gene_dictionary, orient="index")
#print(genes_df.head())

# Load genome dictionary and change to dataframe
with open("/storage/enyaa/REVISED/KMER/genome_dist/genome_kmer_distributions_1.pkl", "rb") as file:
    genome_dictionary = pickle.load(file)

genomes_df = pd.DataFrame.from_dict(genome_dictionary, orient="index")

distances = []
# Loop through matching gene-genome pairs
for _, row in not_mobile_taxonomy.iterrows():
    gene_name = row["Gene_name"]
    bacteria_id = row["Bacteria_ID"]

    # Ensure both k-mer distributions exist
    if gene_name in genes_df.index and bacteria_id in genomes_df.index:
        gene_kmers = genes_df.loc[gene_name]
        genome_kmers = genomes_df.loc[bacteria_id]

        # Compute Euclidean distance
        dist = euclidean(gene_kmers, genome_kmers)
        distances.append({"Gene_name": gene_name, "Bacteria_ID": bacteria_id, "Euclidean_distance": dist})
    

# Convert results to DataFrame
distance_df = pd.DataFrame(distances)
#print(distance_df.head())


plt.figure(figsize=(10, 6))
sns.histplot(distance_df["Euclidean_distance"], bins=30, kde=True, color="blue")

# Labels and title
plt.xlabel("Euclidean distance")
plt.ylabel("Number of bacteria")
plt.title("")

plt.savefig("/home/jolunds/newtest/compatible_references.png")
plt.close()


comp_euclidean_df = []  # Initialize an empty list to store results
for _, row in compatible_df.iterrows():
    gene = row['Gene_name']
    bacteria_id = row['Bacteria_ID']
    
    # Get the Euclidean distance from euclidean_df
    if gene in euclidean_df.index and bacteria_id in euclidean_df.columns:
        euclidean_distance = euclidean_df.at[gene, bacteria_id]
    else:
        euclidean_distance = None  # In case there's no matching entry
    
    # Append the results (gene, bacteria_id, and Euclidean distance) to the list
    comp_euclidean_df.append([gene, bacteria_id, euclidean_distance])

# Convert the list to a DataFrame for better readability
comp_euclidean_df = pd.DataFrame(comp_euclidean_df, columns=['Gene_name', 'Bacteria_ID', 'Euclidean_distance'])

# Convert the filtered rows to a DataFrame
temp_incomp_euclidean_df = pd.DataFrame(incomp_euclidean)

valid_columns = [id for id in incomp_bacteria_id if id in temp_incomp_euclidean_df.columns]

# Filter so only 80 genomes left
filtered_df = temp_incomp_euclidean_df[valid_columns]
random_bacteria_ids = np.random.choice(filtered_df.columns, size=100, replace=False)
filtered_100_df = filtered_df[random_bacteria_ids]

print(len(filtered_100_df.columns))
# Filter so we only take 150 genes
random.seed(42)
random_incomp_genes = random.sample(list(temp_incomp_euclidean_df.index), 150)

incomp_euclidean_df = filtered_100_df.loc[random_incomp_genes]
incomp_euclidean_df = incomp_euclidean_df.reset_index().melt(id_vars=["index"], var_name="Bacteria_ID", value_name="Euclidean_distance")
incomp_euclidean_df.rename(columns={'index': 'Gene_name'}, inplace=True)


print(len(incomp_euclidean_df))

'''