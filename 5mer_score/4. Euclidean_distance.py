# Create dictionary for euclidean distance (5mer score) between genes and genomes

import pandas as pd
import pickle
from scipy.spatial.distance import cdist
import os
from tqdm import tqdm


# Gene kmer distributions
with open ("/storage/enyaa/FINAL/KMER/gene_kmer_distributions.pkl", "rb") as file:
    gene_dictionary = pickle.load(file)
genes_df = pd.DataFrame.from_dict(gene_dictionary, orient="index").T
    # 6048 genes

# Genome kmer distributions
with open ("/storage/enyaa/FINAL/KMER/genome_kmer_distributions.pkl", "rb") as file_genomes:
    genome_dictionary = pickle.load(file_genomes)
genomes_df = pd.DataFrame.from_dict(genome_dictionary, orient="index").T.fillna(0)
    # 72690 genomes

# Fix so both have same index
genes_df = genes_df.reindex(genomes_df.index).fillna(0) 
    
# Calculate the euclidean distance between the two distributions
euclidean = cdist(genes_df.T, genomes_df.T, metric='euclidean')
    
euclidean_df = pd.DataFrame(euclidean, 
                           index=genes_df.columns, 
                           columns=genomes_df.columns)
    

# Save each row (gene) as a separate .pkl file
output_dir = "/storage/enyaa/FINAL/KMER/euclidean_split_genes/"
for gene_name in tqdm(euclidean_df.index, desc="Processing genes"):
    gene_df = euclidean_df.loc[[gene_name]] 
    
    gene_df = gene_df.reset_index(drop=True).T.reset_index() # switch to long format
    gene_df.columns = ['Bacteria_ID', 'Euclidean_distance']
   
    # Add taxonomic info and match status
    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    info_df = pd.read_csv(f"/storage/jolunds/FINAL/MATCH_INFO/{gene_name}_match_info.tsv", sep="\t")
        # 72690 rows (bacteria), columns: Bacteria_ID, Domain, Phylum, Class, Order, Family, Genus, Species, Match_info

    gene_df = pd.merge(gene_df, info_df, on='Bacteria_ID', how='left')    

    # Save
    gene_df.to_pickle(os.path.join(output_dir, f"euclidean_df_{gene_name}.pkl"))

print("Done! Each gene has its own .pkl file.")

