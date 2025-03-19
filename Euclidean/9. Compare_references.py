import json
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import euclidean
import random
import numpy as np

# COMPATIBLE -------------------------------------
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
        
#print(len(chromosomal_genes)) #995

# Load mobility
path = "/storage/enyaa/REVISED/mobility_classification_all.csv"
mobility_df = pd.read_csv(path, sep=",")

# Filter for the chromosomal genes
chromosomal_df = mobility_df[mobility_df['Gene_name'].isin(chromosomal_genes)]
#print(len(chromosomal_df))

# Keep only the "not mobile" genes
not_mobile_df = chromosomal_df[chromosomal_df['Mobility'] == 'Not_mobile']

# Load taxonomy results 1 & 2
file_path_1 = "/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_1.csv"
taxonomy_df_1 = pd.read_csv(file_path_1, sep=",", header=None) # create a pandas dataframe
taxonomy_df_1.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
file_path_2 = "/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_2.csv"
taxonomy_df_2 = pd.read_csv(file_path_2, sep=",", header=None) # create a pandas dataframe
taxonomy_df_2.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

taxonomy_df = pd.concat([taxonomy_df_1, taxonomy_df_2], ignore_index=True) # Merge

not_mobile_taxonomy = not_mobile_df[['Gene_name']].merge(taxonomy_df, on='Gene_name', how='left') #Take 

compatible_df = not_mobile_taxonomy[['Gene_name', 'Bacteria_ID']] # Keep only gene name and bacteria id

# Load euclidean df
#euclidean_df = pd.read_pickle("/storage/enyaa/REVISED/KMER/euclidean_df.pkl") # NOT CREATED YET
euclidean_df = pd.read_pickle("/storage/jolunds/euclidean_df.pkl") 

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
comp_euclidean_df['Reference'] = 'Compatible' # Add 'Compatible'

# INCOMPATIBLE ----------
gram_positive_genera = ["Clostridium", "Bacillus", "Clostridioides", "Listeria", "Staphylococcus", "Enterococcus", "Lactobacillus", "Leuconostoc",
                        "Streptococcus"]

# Load full taxonomy and filter for gram_positive genera
path = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_lineage_df = pd.read_csv(path, sep="\t", header=None)
full_lineage_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

gram_positive_bacteria = full_lineage_df[full_lineage_df["Genus"].isin(gram_positive_genera)]
incomp_bacteria_id = gram_positive_bacteria["Bacteria_ID"]

incomp_gene_names = ["NDM", "IMP", "GIM", "SPM", "VIM"] #Metallo betalaktamases 

incomp_euclidean = []
for idx, row in euclidean_df.iterrows():
    gene = str(idx).strip()  # Ensure it's a string and strip whitespace
    if any(gene.startswith(prefix) for prefix in incomp_gene_names):
        # If the condition is met, append the row to incomp_euclidean
        incomp_euclidean.append(row)
 
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
incomp_euclidean_df['Reference'] = 'Incompatible'

print(len(incomp_euclidean_df))

# Combine incomp and comp
reference_euclidean_df = pd.concat([comp_euclidean_df, incomp_euclidean_df], ignore_index=True)

# Plot results
plt.figure(figsize=(10, 6))
sns.histplot(data=reference_euclidean_df, x='Euclidean_distance', hue='Reference', multiple='stack', bins=20)

plt.xlabel("Euclidean distance")
plt.ylabel("Number of bacteria")

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

'''