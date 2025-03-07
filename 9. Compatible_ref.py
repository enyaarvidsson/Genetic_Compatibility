# Super-compatible reference gene


# !!!!!!!!!!
# BEHÖVER KOLLA IGENOM PGA GENE_NAME OCH ANNAT
# !!!!!!!!!!


# Gene ID: AE006468.2 - golS - chromosomal, in bacteria "Salmonella enterica"
# NG_050391.1 - CfiA2 - in "Bacteroides fragilis"

# ARO:3002571 - chromosomal see CARD website
# ARO:0010002 - chromosomal see reference in google docs

import pickle
import random 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

'''
# för att hitta gene_name från gene_id

def find_gene_name(gene_id, file_path="/storage/enyaa/gene_ids_and_names.txt"):
    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")  # Split by tab
            if parts[0] == gene_id:  # Compare with gene_id
                return parts[1]  # Return gene_name
    return "Gene ID not found"

# Example usage
gene_id_to_search = "NG_050391.1"  # Replace with the actual Gene ID you're looking for
gene_name = find_gene_name(gene_id_to_search)
print(f"Gene Name for {gene_id_to_search}: {gene_name}")

'''
# Read gene dictionary
with open("/storage/enyaa/REVISED/KMER/gene_dist/gene_kmer_distributions.pkl", "rb") as file: #"rb": read binary
    gene_dictionary = pickle.load(file)

# Find the desired gene
gene_id = "NG_050391.1" #  Gene ID for a gene that we picked (chromosomal)
gene_distribution = gene_dictionary.get(gene_id) # Get the desired gene from the dictionary

if gene_distribution is not None:
    print(f"Gene {gene_id} found in the dictionary.")
else:
    print(f"Gene {gene_id} not found in the dictionary.")

# Read genome dictionary
with open("/storage/enyaa/REVISED/KMER/genome10k_kmer_distributions.pkl", "rb") as file:
    genome_dictionary = pickle.load(file)

# Take out genomes that matches with gene AE006468.2
gene_name = "CfiA2"
path_taxonomy = "/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_4.csv"
taxonomy_df = pd.read_csv(path_taxonomy, sep=",", header=None)
filtered_df = taxonomy_df[taxonomy_df.iloc[:,0] == gene_name]   # Takes out the information for specific gene
filtered_df.columns = ["Gene_ID", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
print(filtered_df)

chromosomal_bacteria = ["Bacteroides fragilis"] 
chromosomal_bacteria_df = filtered_df[filtered_df["Species"].isin(chromosomal_bacteria)]

bacteria_ids = chromosomal_bacteria_df["Bacteria_ID"].unique().tolist()
print(bacteria_ids)  

genomes_distributions = {bacteria_id: genome_dictionary[bacteria_id] for bacteria_id in bacteria_ids} # Uses random_bacteria_ids to find 
    # the right keys in the dictionary and takes out the distributions for them 

# Make gene and genome distrbutions into dataframe
gene_df = pd.DataFrame(list(gene_distribution.items()), columns=["kmer", f"{gene_id}"])
gene_df.set_index("kmer", inplace=True)
genomes_df = pd.DataFrame.from_dict(genomes_distributions, orient="index").T

# Join gene with genomes
distributions_df = gene_df.join(genomes_df, how="outer")
distributions_df.fillna(0, inplace=True)

# Calculate eucldian distances between gene and genomes distributions
gene_vector = distributions_df[f"{gene_id}"].values[:, None]
eu_distances = np.linalg.norm(distributions_df.iloc[:, 1:].values - gene_vector, axis=0)

# Store bacteria id with the euclidiean distances 
euclidean_df = pd.DataFrame({
    'Bacteria_ID': distributions_df.columns[1:],  # The genome names (exclude the first column)
    'Euclidean_distance': eu_distances
})
euclidean_df["Gene_ID"] = gene_id

save_path = "/storage/enyaa/compatible_reference.tsv"
euclidean_df.to_csv(save_path, sep="\t", index=False)



# Create histogram with stacked bars
plt.figure(figsize=(8, 5))

#palette = sns.color_palette(["#1F77B4", "#FF7F0E"])  # orange for "Match" and blue for "No Match"
# Apply the palette to the plot
#sns.set_palette(palette)

sns.histplot(data=euclidean_df, x='Euclidean_distance', multiple='stack',  bins=5) 

plt.xlabel('Euclidean distance')
plt.ylabel('Number of genomes')
plt.title('Title')

plt.savefig('/home/enyaa/gene_genome/histogram_compatible.png') 
plt.close()

