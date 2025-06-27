import pandas as pd
import pickle
import numpy as np
import time
import os
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import seaborn as sns

start_time = time.time()

# Load k-mer distributions for all genes    
with open ("/storage/enyaa/FINAL/KMER/gene_kmer_distributions.pkl", "rb") as file:
    gene_dictionary = pickle.load(file)

genes_df = pd.DataFrame.from_dict(gene_dictionary, orient="index").T #Change to dataframe

file_path = f"/storage/enyaa/FINAL/KMER/genome_kmer_distributions.pkl"
with open (file_path, "rb") as file:
    genome_dictionary = pickle.load(file)
    
genomes_df = pd.DataFrame.from_dict(genome_dictionary, orient="index").T.fillna(0) # Change to dataframe
        
genes_df = genes_df.reindex(genomes_df.index).fillna(0) # Fill missing k-mers in genes_df

# Compute max difference row by row
for gene_name in genes_df.columns:
    gene_values = genes_df[gene_name].values  # Extract 1D array 
        
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
            
    path = f"/storage/jolunds/FINAL/WORST_CASE/worst_case_{gene_name}.csv"
    worst_case_df.to_csv(path)



# Scatterplot 5mer score vs 5mer worst case combined
gene_name = "tet(Q)" 

np.random.seed(42)

if "/" in gene_name:
    gene_name = gene_name.replace("/", "?")


# EUCLIDEAN DISTANCE ------------------
filepath_eu = f"/storage/enyaa/FINAL/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl"
euclidean_df = pd.read_pickle(filepath_eu)

no_match_count = euclidean_df["Match_status"].value_counts().get("Match", 0)    
if no_match_count == 0:
    #print("No matches for gene:", gene_name)
    matches = 0

# DOWNSAMPLE NO-MATCHES -------
euclidean_downsampled = [] # will become a list of dataframes

match_count = (euclidean_df["Match_status"] == "Match").sum()

# create a df for match and a df for no_match
matches_df = euclidean_df[euclidean_df["Match_status"] == "Match"]
no_matches_df = euclidean_df[euclidean_df["Match_status"] == "Non_match"]  

# How many no_matches to keep
if matches == 1: # if matches exists
    if match_count < 100:
        keep_size = 100
    else:
        keep_size = match_count
else:
    print("No matches") 

# Downsample no_matches
if keep_size > len(no_matches_df): 
    keep_size = len(no_matches_df)
downsampled_no_matches_df = no_matches_df.iloc[np.random.choice(len(no_matches_df), keep_size, replace=False)]

# Append both "Match" and downsampled "Non_match" bacteria
euclidean_downsampled.append(pd.concat([matches_df, downsampled_no_matches_df]))

euclidean_downsampled_df = pd.concat(euclidean_downsampled, ignore_index=True)

# WORST CASE ------------------
file_path = f"/storage/jolunds/FINAL/WORST_CASE/worst_case_{gene_name}.csv"
worst_case_df = pd.read_csv(file_path)
worst_case_df = worst_case_df.drop(worst_case_df.columns[0], axis=1)

# MERGE THEM TOGETHER ----------
euclidean_and_worst_df = pd.merge(euclidean_downsampled_df, worst_case_df, on='Bacteria_ID', how='inner')


# WEIGHTED COMBINATION OF MAX DIFF AND RELATIVE DIFF -----------
scaler = MinMaxScaler()
euclidean_and_worst_df[["Max_diff_scaled", "Rel_diff_scaled"]] = scaler.fit_transform(
    euclidean_and_worst_df[["Max_difference", "Relative_difference"]]
)

euclidean_and_worst_df["Combined_score"] = (
    0.5 * euclidean_and_worst_df["Max_diff_scaled"] +
    0.5 * euclidean_and_worst_df["Rel_diff_scaled"]
)

# SCATTERPLOT -----------
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=euclidean_and_worst_df,
    x='Combined_score',
    y='Euclidean_distance',   
    hue='Match_status',
    hue_order=["Non-match", "Match"],
    alpha=1,
    s=20
    #color='darkorange'
)

# Add title
if "?" in gene_name:
    gene_name = gene_name.replace("?", "/")
 
plt.xlabel('5mer worst case', fontsize=14)
plt.ylabel('5mer score', fontsize=14)
plt.tick_params(axis='both', labelsize=12)
plt.legend(fontsize=14, loc="upper center")
plt.tight_layout()

plt.savefig(f'./5mer_vs_worst_{gene_name}.png')     
plt.close()

print(f"Scatterplot created for {gene_name}")



