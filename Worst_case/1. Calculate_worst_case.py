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


max_diff_data = []
rel_diff_data = []

for i in range(1, 17):
    file_path = f"/storage/enyaa/REVISED/KMER/genome_dist/genome_kmer_distributions_{i}.pkl"
    with open(file_path, "rb") as file:
        genome_dictionary = pickle.load(file)

    genomes_df = pd.DataFrame.from_dict(genome_dictionary, orient="index").T.fillna(0)
    genes_df = genes_df.reindex(genomes_df.index).fillna(0)

    for gene in genes_df.columns:
        gene_values = genes_df[gene].values

        abs_diff = np.abs(gene_values[:, np.newaxis] - genomes_df.values)
        max_diff = np.max(abs_diff, axis=0)

        mean_values = (gene_values[:, np.newaxis] + genomes_df.values) / 2
        rel_diff = np.max(abs_diff / mean_values, axis=0)

        # Append results directly to lists
        max_diff_data.extend(zip([gene] * len(genomes_df.columns), genomes_df.columns, max_diff))
        rel_diff_data.extend(zip([gene] * len(genomes_df.columns), genomes_df.columns, rel_diff))

    print(f"Done with database {i}")

# Convert lists to DataFrames (faster than continuous DataFrame appends)
max_diff_df = pd.DataFrame(max_diff_data, columns=["Gene_name", "Bacteria_ID", "Max_Diff"])
rel_diff_df = pd.DataFrame(rel_diff_data, columns=["Gene_name", "Bacteria_ID", "Rel_Diff"])

# Merge them efficiently
worst_case_df = pd.concat([max_diff_df, rel_diff_df["Rel_Diff"]], axis=1)
print(worst_case_df.head())

output_directory = f"/storage/jolunds/REVISED/WORST_CASE/new_worst_case_split/"

for gene_name, gene_df in worst_case_df.groupby("Gene_name"):
    file_path = os.path.join(output_directory, f"worst_case_{gene_name}.pkl")
    gene_df.to_pickle(file_path)


'''num_files = 0
output_directory = f"/storage/jolunds/REVISED/WORST_CASE/worst_case_split/"

for row in worst_case_df.itertuples(index=True, name=None):
    gene_name = row[0]
    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")
        
    filename = os.path.join(output_directory, f"worst_case_{gene_name}.csv")
    # Save the row as a single-row CSV
    with open(filename, "w") as f:
        f.write(",".join(worst_case_df.columns) + "\n")  # Write header
        f.write(",".join(map(str, row[1:])) + "\n")  # Write values
    
    num_files += 1
'''
end_time = time.time()
total_time = (end_time - start_time)/60 

#print(f"Done creating {num_files} files in {total_time} minutes")

