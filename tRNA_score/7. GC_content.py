# This script includes three parts:
    # Calculating the GC-ratio and making a scatterplot (GC-ratio vs tRNA score)
    # Calculating the GC-difference and making a scatterplot (GC-difference vs tRNA score)
    # For the GC-ratio plot, dividing it up into bins and making all bins have the same number of datapoints

import time
import pandas as pd
import pickle
import os
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt


# SCATTERPLOT RATIO for matches -----------------------------
# Calculates the GC-ratio between genes and genomes - only for filtered matches
# Makes a scatterplot gc-ratio vs tRNA score

""" start_time = time.time()

# Load the files with GC-content for genes and genomes 
file_bacteria = "/storage/enyaa/FINAL/GC/gc_content_bacteria.pkl"
file_genes = "/storage/enyaa/FINAL/GC/gc_content_genes.pkl"

with open(file_bacteria, "rb") as f:
    bacteria_gc_df = pickle.load(f)
with open(file_genes, "rb") as f:
    genes_gc_df = pickle.load(f)


# Go through each gene
gene_names = "/storage/enyaa/FINAL/gene_names.txt"
gene_names_df = pd.read_csv(gene_names, header=None, names=["Gene_name"])

tRNA_score_all = []
gc_ratio_all = []

for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): 

    gene_gc = genes_gc_df.loc[genes_gc_df['Gene_name'] == gene_name, 'GC_content'].values[0]
    
    # ONLY COMPUTE GC-RATIO BETWEEN GENES AND ITS MATCHES

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    tRNA_gene_df = pd.read_csv(f"/storage/jolunds/FINAL/tRNA_SCORE/{gene_name}_tRNA_score.csv") 

    # filter to only include matching bacteria
    filtered_tRNA_gene_df = tRNA_gene_df[tRNA_gene_df["Match_status"] == "Match"]

    if filtered_tRNA_gene_df.empty: 
        continue
    
    filtered_gc_df = bacteria_gc_df.merge(
        filtered_tRNA_gene_df[['Bacteria_ID']], on='Bacteria_ID', how='inner'
    )

    bacteria_gc = filtered_gc_df['GC_content'].to_numpy()

    # Compute GC ratio
    ratio = np.round(gene_gc / bacteria_gc, 4) # array of the ratio for one gene vs all genomes

    tRNA_score_all.extend(filtered_tRNA_gene_df["tRNA_score"].values) 
    gc_ratio_all.extend(ratio)

df_plot = pd.DataFrame({
    'tRNA_score': tRNA_score_all,
    'GC_ratio': gc_ratio_all
})

# Scatterplot:
plt.figure(figsize=(8, 6))
plt.scatter(data=df_plot, x='tRNA_score', y='GC_ratio', alpha=1, s=10, color= "darkorange")
plt.xlabel("tRNA score", fontsize=14)
plt.ylabel("GC-ratio", fontsize=14)
plt.tick_params(axis='both', labelsize=12)

#plt.title(f"GC-ratio vs tRNA score for matching genes and genomes")
plt.tight_layout()
plt.grid(True)
plt.savefig(f'/home/enyaa/gene_genome/GC_ratio_vs_tRNA.png') 
plt.close()

end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot ratio created in: {total_time} minutes!") """


# SCATTERPLOT DIFFERENCE for matches -------------------------
# Calculates the GC-difference between genes and genomes - only for filtered matches
# Makes a scatterplot gc-difference vs tRNA score
""" 
start_time = time.time()

# Load the files with GC-content for genes and genomes 
file_bacteria = "/storage/enyaa/FINAL/GC/gc_content_bacteria.pkl"
file_genes = "/storage/enyaa/FINAL/GC/gc_content_genes.pkl"

with open(file_bacteria, "rb") as f:
    bacteria_gc_df = pickle.load(f)
with open(file_genes, "rb") as f:
    genes_gc_df = pickle.load(f)


# Go through each gene
gene_names = "/storage/enyaa/FINAL/gene_names.txt"
gene_names_df = pd.read_csv(gene_names, header=None, names=["Gene_name"])

tRNA_score_all = []
gc_diff_all = []

for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): 

    gene_gc = genes_gc_df.loc[genes_gc_df['Gene_name'] == gene_name, 'GC_content'].values[0]
    
    # ONLY COMPUTE GC-DIFFERENCE BETWEEN GENES AND ITS MATCHES

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    tRNA_gene_df = pd.read_csv(f"/storage/jolunds/FINAL/tRNA_SCORE/{gene_name}_tRNA_score.csv") 

    # filter to only include matching bacteria
    filtered_tRNA_gene_df = tRNA_gene_df[tRNA_gene_df["Match_status"] == "Match"]

    if filtered_tRNA_gene_df.empty: 
        continue
    
    filtered_gc_df = bacteria_gc_df.merge(
        filtered_tRNA_gene_df[['Bacteria_ID']], on='Bacteria_ID', how='inner'
    )

    bacteria_gc = filtered_gc_df['GC_content'].to_numpy()

    # Compute GC difference
    diff = np.round(gene_gc - bacteria_gc, 4) # array of the difference for one gene vs all genomes

    tRNA_score_all.extend(filtered_tRNA_gene_df["tRNA_score"].values) 
    gc_diff_all.extend(diff)

df_plot = pd.DataFrame({
    'tRNA_score': tRNA_score_all,
    'GC_diff': gc_diff_all
})

# Scatterplot:
plt.figure(figsize=(8, 6))
plt.scatter(data=df_plot, x='tRNA_score', y='GC_diff', alpha=1, s=10, color= "darkorange")
plt.xlabel("tRNA score", fontsize=14)
plt.ylabel("GC-difference", fontsize=14)
plt.tick_params(axis='both', labelsize=12)

#plt.title(f"GC-difference vs tRNA score for matching genes and genomes")
plt.tight_layout()
plt.grid(True)
plt.savefig(f'/home/enyaa/gene_genome/GC_diff_vs_tRNA.png') 
plt.close()

end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot difference created in: {total_time} minutes!") """


# SCATTERPLOT RATIO - BINS for matches -----------------------------
# Calculates the GC-ratio between genes and genomes - only for filtered matches
# Divides into bins
# Makes a scatterplot gc-ratio vs tRNA score


start_time = time.time()

# Load the files with GC-content for genes and genomes 
file_bacteria = "/storage/enyaa/FINAL/GC/gc_content_bacteria.pkl"
file_genes = "/storage/enyaa/FINAL/GC/gc_content_genes.pkl"

with open(file_bacteria, "rb") as f:
    bacteria_gc_df = pickle.load(f)
with open(file_genes, "rb") as f:
    genes_gc_df = pickle.load(f)

# Go through each gene
gene_names = "/storage/enyaa/FINAL/gene_names.txt"
gene_names_df = pd.read_csv(gene_names, header=None, names=["Gene_name"])

tRNA_score_all = []
gc_ratio_all = []

for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): 

    gene_gc = genes_gc_df.loc[genes_gc_df['Gene_name'] == gene_name, 'GC_content'].values[0]
    
    # ONLY COMPUTE GC-RATIO BETWEEN GENES AND ITS MATCHES

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    tRNA_gene_df = pd.read_csv(f"/storage/jolunds/FINAL/tRNA_SCORE/{gene_name}_tRNA_score.csv") 

    # filter to only include matching bacteria
    filtered_tRNA_gene_df = tRNA_gene_df[tRNA_gene_df["Match_status"] == "Match"]

    if filtered_tRNA_gene_df.empty: 
        continue
    
    filtered_gc_df = bacteria_gc_df.merge(
        filtered_tRNA_gene_df[['Bacteria_ID']], on='Bacteria_ID', how='inner'
    )

    bacteria_gc = filtered_gc_df['GC_content'].to_numpy()

    # Compute GC ratio
    ratio = np.round(gene_gc / bacteria_gc, 4) # array of the ratio for one gene vs all genomes

    tRNA_score_all.extend(filtered_tRNA_gene_df["tRNA_score"].values) 
    gc_ratio_all.extend(ratio)

df_plot = pd.DataFrame({
    'tRNA_score': tRNA_score_all,
    'GC_ratio': gc_ratio_all
})


# BINS ---
# Add a column bin to the df, with the labels
bins = [0, 0.5, 0.75, 1.0, 1.25, 1.5, float('inf')]
labels = ['<0.5', '0.5-0.75', '0.75-1.0', '1.0-1.25', '1.25-1.5', '>1.5']
df_plot['bin'] = pd.cut(df_plot['GC_ratio'], bins=bins, labels=labels, right=False)

# Only want to look at the following bins at the moment
target_bins = ['0.5-0.75', '0.75-1.0', '1.0-1.25', '1.25-1.5']

# Filter the df to only include rows from the target bins
df_target = df_plot[df_plot['bin'].isin(target_bins)]

# Create a copy to avoid the SettingWithCopyWarning
df_target = df_target.copy()

# Drop labels not used (<0.5 and >1.5)
df_target['bin'] = df_target['bin'].cat.remove_unused_categories()

# Count how many in each bin
bin_counts = df_target['bin'].value_counts()
print(bin_counts)

# Get the smallest count
min_count = bin_counts.min()
print(min_count)

# Group by bin and sample min_count rows from each group
sampled_dfs = (
    df_target
    .groupby('bin', group_keys=False, observed=True)
    .apply(lambda x: x.sample(n=min_count, random_state=42))
    .reset_index(drop=True)
)

# Keep all rows from <0.5 and >1.5 
other_bins_df = df_plot[df_plot['bin'].isin(['<0.5', '>1.5'])]

# Concatenate the result
balanced_df = pd.concat([sampled_dfs, other_bins_df], ignore_index=True)

df_plot = balanced_df

# Scatterplot:
plt.figure(figsize=(8, 6))
plt.scatter(data=df_plot, x='tRNA_score', y='GC_ratio', alpha=1, s=10, color= "darkorange")
plt.xlabel("tRNA score", fontsize=14)
plt.ylabel("GC-ratio", fontsize=14)
plt.tick_params(axis='both', labelsize=12)


#plt.title(f"GC-ratio vs tRNA score for matching genes and genomes")
#plt.legend(title='Phylum', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.grid(True)
plt.savefig(f'/home/enyaa/gene_genome/bins.png') 
plt.close()


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot ratio bins created in: {total_time} minutes!")

