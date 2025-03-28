import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import time
import os

start_time = time.time()
    
with open ("/storage/enyaa/REVISED/KMER/gene_dist/gene_kmer_distributions.pkl", "rb") as file:
    gene_dictionary = pickle.load(file)

#gene_dictionary = dict(list(gene_dictionary.items())[:1])
genes_df = pd.DataFrame.from_dict(gene_dictionary, orient="index").T

for i in range(1,17):
    file_path = f"/storage/enyaa/REVISED/KMER/genome_dist/genome_kmer_distributions_{i}.pkl"
    with open (file_path, "rb") as file:
        genome_dictionary = pickle.load(file)
    
    #genome_dictionary = dict(list(genome_dictionary.items())[:1000])
    genomes_df = pd.DataFrame.from_dict(genome_dictionary, orient="index").T.fillna(0)
        
    genes_df = genes_df.reindex(genomes_df.index).fillna(0) 

    # Compute max difference row by row
    for gene in genes_df.columns:
        gene_values = genes_df[gene].values  # Extract 1D array (size 512)
        # Compute absolute differences
        abs_diff = np.abs(gene_values[:, np.newaxis] - genomes_df.values)

        # Find max difference for each column
        max_diff = np.max(abs_diff, axis=0)

        # Find indices where max difference occurs
        max_indices = np.argmax(abs_diff, axis=0)

        # Extract values that gave max difference
        gene_max_values = gene_values[max_indices]
        #print(gene_max_values)
        genome_max_values = genomes_df.values[max_indices, np.arange(genomes_df.shape[1])]
        #print(genome_max_values)
        # Compute relative difference
        relative_diff = max_diff / ((np.abs(gene_max_values + genome_max_values)) / 2)
        
        # Store results in dataframe
        worst_case_df = pd.DataFrame({
        "Bacteria_ID": genomes_df.columns,  # Assuming columns are bacteria names
        "Max_difference": max_diff,
        "Relative_difference": relative_diff
        })
        
        #print(worst_case_df.head())
        # Save file
        if "/" in gene:
            gene = gene.replace("/", "?")
            
        path = f"/storage/jolunds/REVISED/WORST_CASE/new_worst_case_split/worst_case_{i}_{gene}.csv"
        worst_case_df.to_csv(path)
        

    
    print(f"Done with database {i}")
