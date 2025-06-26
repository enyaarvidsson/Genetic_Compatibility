# Create dictionary for euclidean distance (5mer score) between 
# genomes and length-adjusted genes 

import pandas as pd
import pickle
from scipy.spatial.distance import cdist
from tqdm import tqdm


# Gene kmer distributions
with open ("/storage/enyaa/FINAL/KMER/FOR_GENE_LENGTH/gene_kmer_distributions_500bp.pkl", "rb") as file:
    gene_dictionary = pickle.load(file)
genes_df = pd.DataFrame.from_dict(gene_dictionary, orient="index").T
    # 5887 genes

# Genome kmer distributions
with open ("/storage/enyaa/FINAL/KMER/genome_kmer_distributions.pkl", "rb") as file_genomes:
    genome_dictionary = pickle.load(file_genomes)
genomes_df = pd.DataFrame.from_dict(genome_dictionary, orient="index").T.fillna(0)
    # 72690 genomes

# Fix so both have same index
genes_df = genes_df.reindex(genomes_df.index).fillna(0) 
    
# Calculated the euclidean distance between the two distributions
euclidean = cdist(genes_df.T, genomes_df.T, metric='euclidean')
    
euclidean_df = pd.DataFrame(euclidean, 
                           index=genes_df.columns, 
                           columns=genomes_df.columns)

# Save 
euclidean_df.to_pickle("/storage/enyaa/FINAL/KMER/FOR_GENE_LENGTH/euclidean_df_500bp.pkl")

print("Done calculating euclidean distance for length-adjusted genes!")

