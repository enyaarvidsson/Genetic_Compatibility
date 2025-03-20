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

genes_df = pd.DataFrame.from_dict(gene_dictionary, orient="index").T

worst_case_list = []    
for i in range(1,3):
    file_path = f"/storage/enyaa/REVISED/KMER/genome_dist/genome_kmer_distributions_{i}.pkl"
    with open (file_path, "rb") as file:
        genome_dictionary = pickle.load(file)
    
    genomes_df = pd.DataFrame.from_dict(genome_dictionary, orient="index").T.fillna(0)
        
    genes_df = genes_df.reindex(genomes_df.index).fillna(0) 
    
    max_difference_df = pd.DataFrame(index=genes_df.columns, columns=genomes_df.columns)
    
    # Compute max difference row by row
    for gene in genes_df.columns:
        gene_values = genes_df[gene].values  # Extract 1D array (size 512)
    
        # Compute absolute difference row-wise, then take max
        max_difference_df.loc[gene] = np.max(np.abs(gene_values[:, np.newaxis] - genomes_df.values), axis=0)

    # Convert to numeric values
    max_difference_df = max_difference_df.apply(pd.to_numeric)
    
    
    worst_case_list.append(max_difference_df)
    print(f"Done with database {i}")

worst_case_df = pd.concat(worst_case_list, axis=1)

num_files = 0
output_directory = f"/storage/jolunds/REVISED/WORST_CASE/worst_case_split/"
for gene_name, row in worst_case_df.iterrows():
    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")
        
    filename = os.path.join(output_directory, f"worst_case_{gene_name}.csv")
    row.to_frame().T.to_csv(filename, index=False)
    
    num_files += 1

end_time = time.time()
total_time = (end_time - start_time)/60 

print(f"Done creating {num_files} files in {total_time} minutes")

