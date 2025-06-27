import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import time
import os

'''start_time = time.time()

# Load k-mer distributions for all genes    
with open ("/storage/enyaa/REVISED/KMER/gene_dist/gene_kmer_distributions.pkl", "rb") as file:
    gene_dictionary = pickle.load(file)

genes_df = pd.DataFrame.from_dict(gene_dictionary, orient="index").T #Change to dataframe

# Loop through the files of k-mer distributions for the genomes
for i in range(1,17):
    file_path = f"/storage/enyaa/REVISED/KMER/genome_dist/genome_kmer_distributions_{i}.pkl"
    with open (file_path, "rb") as file:
        genome_dictionary = pickle.load(file)
    
    genomes_df = pd.DataFrame.from_dict(genome_dictionary, orient="index").T.fillna(0) # Change to dataframe
        
    genes_df = genes_df.reindex(genomes_df.index).fillna(0) # Fill missing k-mers in genes_df

    # Compute max difference row by row
    for gene in genes_df.columns:
        gene_values = genes_df[gene].values  # Extract 1D array 
        
        # Compute absolute differences
        abs_diff = np.abs(gene_values[:, np.newaxis] - genomes_df.values)

        # Find max difference for each column
        max_diff = np.max(abs_diff, axis=0)

        # Find indices where max difference occurs
        max_indices = np.argmax(abs_diff, axis=0)

        # Extract values that gave max difference
        gene_max_values = gene_values[max_indices]
        genome_max_values = genomes_df.values[max_indices, np.arange(genomes_df.shape[1])]

        # Compute relative difference
        relative_diff = max_diff / ((np.abs(gene_max_values + genome_max_values)) / 2)
        
        # Store results in dataframe
        worst_case_df = pd.DataFrame({
        "Bacteria_ID": genomes_df.columns,  
        "Max_difference": max_diff,
        "Relative_difference": relative_diff
        })
        
        # Save file
        if "/" in gene:
            gene = gene.replace("/", "?")
            
        path = f"/storage/jolunds/REVISED/WORST_CASE/new_worst_case_split/worst_case_{i}_{gene}.csv"
        worst_case_df.to_csv(path)
        
    print(f"Done with database {i}")
'''
# Merge files
start_time = time.time()
import os
import pandas as pd
from glob import glob

# Paths
directory = "/storage/jolunds/REVISED/WORST_CASE/new_worst_case_split/"
output_dir = "/storage/jolunds/REVISED/WORST_CASE/worst_case_split/"

# Load the list of gene names 
gene_list_file = "/storage/jolunds/REVISED/gene_names.txt"  
with open(gene_list_file, "r") as f:
    gene_names = [line.strip() for line in f]

# Loop through each gene name
for gene_name in gene_names:
    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")
    
    matching_files = glob(os.path.join(directory, f"worst_case_*_{gene_name}.csv")) # Find matching files
    
    if not matching_files:
        print(f"No files found for {gene_name}, skipping.")
        continue

    # Read and merge all CSVs for the current gene
    df_list = [pd.read_csv(file) for file in matching_files]
    merged_df = pd.concat(df_list, ignore_index=True)

    # Save merged file
    output_file = os.path.join(output_dir, f"worst_case_{gene_name}.csv")
    merged_df.to_csv(output_file, index=False)

end_time = time.time()
total_time = (end_time-start_time)/60
print(f"Merging complete, took {total_time} minutes")

