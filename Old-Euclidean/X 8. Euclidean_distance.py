# Create dictionary for euclidean distance between genes and genomes

import pandas as pd
import pickle
from scipy.spatial.distance import cdist
from tqdm import tqdm

#with open ("/storage/enyaa/REVISED/KMER/gene_dist/gene_kmer_distributions.pkl", "rb") as file:
with open ("/storage/enyaa/REVISED/KMER/FOR_GENE_LENGTH/gene_kmer_distributions_500bp.pkl", "rb") as file:
    gene_dictionary = pickle.load(file)

genes_df = pd.DataFrame.from_dict(gene_dictionary, orient="index").T

# Load one euclidean df file to get the filtered bacteria ids (77k)
file = "/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_tet(Q).pkl"
filtered_df = pd.read_pickle(file)
bacteria_ids = filtered_df["Bacteria_ID"].tolist()

euclidean_list_df = []
for i in tqdm(range(1,17), desc="Processing"): # start at first nr, end before second nr
    file_path = f"/storage/enyaa/REVISED/KMER/genome_dist/genome_kmer_distributions_{i}.pkl"
    
    with open (file_path, "rb") as file:
        genome_dictionary = pickle.load(file)
    
    genomes_df = pd.DataFrame.from_dict(genome_dictionary, orient="index").T.fillna(0)
    #print(genomes_df.head())

    # filter columns to keep only those in bacteria_ids 
    genomes_df = genomes_df.loc[:, genomes_df.columns.isin(bacteria_ids)]
    #print(filtered_genomes_df.head())
    
    genes_df = genes_df.reindex(genomes_df.index).fillna(0) 
    #print(genes_df.head())
    
    euclidean_loop = cdist(genes_df.T, genomes_df.T, metric='euclidean')
    
    distance_df = pd.DataFrame(euclidean_loop, 
                           index=genes_df.columns, 
                           columns=genomes_df.columns)
    
    #print(distance_df.head())
    euclidean_list_df.append(distance_df)

euclidean_df = pd.concat(euclidean_list_df, axis=1)

#print(euclidean_df.head())

# Save 
euclidean_df.to_pickle("/storage/enyaa/REVISED/KMER/FOR_GENE_LENGTH/euclidean_df_all_500bp.pkl")
