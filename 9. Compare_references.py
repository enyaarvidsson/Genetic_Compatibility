import json
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import euclidean

# COMPATIBLE ----------
with open("/storage/jolunds/card.json", "r") as file:
    card_info = json.load(file)

# Search for chromosom
search_word = "chromosom"

# Find names of chromosomal genes
chromosomal_genes = []
for key, value in card_info.items():
    if search_word in str(value):
        # Check if 'model_name' exists in the value 
        if isinstance(value, dict) and 'model_name' in value:
           chromosomal_genes.append(value['model_name'])
        
#print(len(chromosomal_genes)) #995

# Load mobility
path = "/storage/jolunds/mobility_classification_1.csv"
mobility_df = pd.read_csv(path, sep=",")

# Filter for the chromosomal genes
chromosomal_df = mobility_df[mobility_df['Gene_name'].isin(chromosomal_genes)]
#print(len(chromosomal_df))

# Keep only the "not mobile" genes
not_mobile_df = chromosomal_df[chromosomal_df['Mobility'] == 'Not_mobile']

# Load taxonomy results 
taxonomy_df = pd.read_csv("/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_1.csv", header=None) # create a pandas dataframe
taxonomy_df.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

not_mobile_taxonomy = taxonomy_df[taxonomy_df['Gene_name'].isin(not_mobile_df['Gene_name'])]
compatible_df= not_mobile_taxonomy[['Gene_name', 'Bacteria_ID']] # Keep only gene name and bacteria id

# INCOMPATIBLE ----------
gram_positive_genera = ["Bacillus", "Listeria", "Staphylococcus", "Enterococcus", "Lactobacillus", "Leuconostoc",
                        "Streptococcus", "Clostridioides", "Clostridium", "Erysipelothrix"]

# Load full taxonomy and filter for gram_positive genera
path = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_lineage_df = pd.read_csv(path, sep="\t", header=None)
full_lineage_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

gram_positive_bacteria = full_lineage_df[full_lineage_df["Genus"].isin(gram_positive_genera)]
incomp_bacteria_id = gram_positive_bacteria["Bacteria_ID"]
incomp_gene_names = ["NDM", "VIM", "IMP", "CphA"] #Metallo betalaktamases 

# BOTH -----------
# Load euclidean df
euclidean_df = pd.read_csv("/storage/enyaa/REVISED/KMER/euclidean_df.pkl") # NOT CREATED YET

comp_euclidean_df = euclidean_df.loc[compatible_df['Gene_name'], compatible_df['Bacteria_ID']] # Compatible referemce

incomp_euclidean_df = euclidean_df.loc[:, euclidean_df.columns.isin(incomp_bacteria_id)] # Incompatible reference
incomp_euclidean_df = incomp_euclidean_df[incomp_euclidean_df.index.str.startswith(tuple(incomp_gene_names))]

# Merge to one 
comp_values = comp_euclidean_df.values.flatten()
incomp_values = incomp_euclidean_df.values.flatten()

df_plot = pd.DataFrame({
    'Euclidean Distance': list(comp_values) + list(incomp_values),
    'Reference': ['Compatible'] * len(comp_values) + ['Incompatible'] * len(incomp_values)
})

# Plot results
plt.figure(figsize=(10, 6))
sns.histplot(data=df_plot, x='Euclidean Distance', hue='Reference', bins=20, kde=True, alpha=0.5)

plt.xlabel("Euclidean distance")
plt.ylabel("Number of bacteria")
plt.title("")

plt.savefig("/home/jolunds/newtest/compatible_references.png")
plt.close()

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