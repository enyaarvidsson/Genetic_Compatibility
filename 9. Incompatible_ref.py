import pandas as pd
from Bio import SeqIO
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Read gene dictionary
with open("/storage/enyaa/REVISED/KMER/gene_dist/gene_kmer_distributions.pkl", "rb") as file:
    gene_dictionary = pickle.load(file)
    
gene_name = "NDM-5" # New Delhi beta-lactamase NDM-5. Chosen from code in the bottom of the script
gene_distribution = gene_dictionary.get(gene_name)

# Find bacteria IDs for gram positive bacteria
path_full_lineage = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_lineage_df = pd.read_csv(path_full_lineage, sep="\t", header=None)
full_lineage_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# List of Gram-positive genera 
gram_positive_genera = ["Clostridium"] #["Bacillus", "Clostridium", "Lactobacillus", "Staphylococcus", "Streptococcus", 
                         #"Enterococcus"]

# Filter for Gram-positive bacteria
gram_positive_bacteria = full_lineage_df[full_lineage_df["Genus"].isin(gram_positive_genera)]
bacteria_id = gram_positive_bacteria["Bacteria_ID"]

# Read genome dictionary
with open("/storage/enyaa/REVISED/KMER/genome_dist/genome_kmer_distributions_1.pkl", "rb") as fil:
    genome_dictionary = pickle.load(fil)

genomes_distributions = {genomes: genome_dictionary.get(genomes) for genomes in bacteria_id} # filter for the bacteria that we have in dictionary
genomes_distributions = {k: v for k, v in genomes_distributions.items() if v is not None} #remove nones

# Calculate euclidean distances
gene_df = pd.DataFrame(list(gene_distribution.items()), columns=["kmer", f"{gene_name}"]) # Change to dataframe
gene_df.set_index("kmer", inplace=True)
genomes_df = pd.DataFrame.from_dict(genomes_distributions, orient="index").T # Change to dataframe

distributions_df = gene_df.join(genomes_df, how="outer") #Merge
distributions_df.fillna(0, inplace=True) # Fill NaNs with 0

gene_vector = distributions_df[f"{gene_name}"].values[:, None]
eu_distances = np.linalg.norm(distributions_df.iloc[:, 1:].values - gene_vector, axis=0)

euclidean_df = pd.DataFrame({
    'Bacteria_ID': distributions_df.columns[1:],  # The genome names (exclude the first column)
    'Euclidean_distance': eu_distances
})
#euclidean_df["Genus"] = euclidean_df["Bacteria_ID"].map(gram_positive_bacteria.set_index("Bacteria_ID")["Genus"]) # Add genus information to euclidean_df
euclidean_df["Gene_name"] = gene_name
#print(euclidean_df.head())

# Double check so there is no matches
file = "/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_1.csv"
taxonomy_results = pd.read_csv(file, sep=",", header=None)
taxonomy_results.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
filtered_df = taxonomy_results[taxonomy_results.iloc[:,0] == gene_name] 

match_column = []
match = 0
no_match = 0
for bacteria in euclidean_df['Bacteria_ID']:
    if bacteria in filtered_df.iloc[:, 1].values:
        match_column.append("Match")
        match += 1
    else:
        match_column.append("No_match")
        no_match += 1
        
print(f"Match: {match}")
print(f"No match: {no_match}")

'''
save_path = "/storage/enyaa/incompatible_reference.tsv"
euclidean_df.to_csv(save_path, sep="\t", index=False)
'''

'''
filepath = "/storage/jolunds/nucleotide_fasta_protein_homolog_model.fasta"
fasta_headers = [record.description for record in SeqIO.parse(filepath, "fasta")]

mbl_genes = ["NDM", "VIM", "IMP", "CphA"]

# Add gene name column in blast_df
classB_betalaktmases = [header for header in fasta_headers if any(mbl_gene in header for mbl_gene in mbl_genes)]

random.seed(42)
random_betalaktamas = random.choice(classB_betalaktmases)

print(random_betalaktamas)
'''

'''
filepath = "/storage/jolunds/shortname_antibiotics.tsv"
shortname_antibiotics_df = pd.read_csv(filepath, sep="\t")

save_path = "/home/jolunds/newtest/shortname_antibiotics.tsv"

shortname_antibiotics_df.to_csv(save_path)

#print(shortname_antibiotics_df)

filepath = "/storage/jolunds/aro_categories_index.tsv"
aro_info = pd.read_csv(filepath, sep="\t")

betalaktamas_df = aro_info[aro_info.apply(lambda row: row.astype(str).str.contains("beta-lactamase", case=False, na=False).any(), axis=1)]
num_rows = len(betalaktamas_df)
print(num_rows)

random_betalaktamas = betalaktamas_df.sample(n=1, random_state=5)

print(random_betalaktamas)
reference_gene = random_betalaktamas.iloc[0,0]

print(reference_gene)

'''