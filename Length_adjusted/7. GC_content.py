# This script includes two parts:
    # Calculating the GC-ratio and making a scatterplot (GC-ratio vs length-adjusted 5mer score)
    # Calculating the GC-difference and making a scatterplot (GC-difference vs length-adjusted 5mer score)

from Bio.SeqUtils import gc_fraction
import pandas as pd
import gzip
from Bio import SeqIO
import os
import pickle
import time
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import seaborn as sns


# SCATTERPLOT RATIO for matches ------------------------------
# Calculates the GC-ratio between genes and genomes - only for filtered matches
# Makes a scatterplot gc-ratio vs length-adjusted 5mer score
 
""" start_time = time.time()

# Load the files with GC-content for genes and genomes
file_bacteria = "/storage/enyaa/FINAL/GC/gc_content_bacteria.pkl"
file_genes = "/storage/enyaa/FINAL/GC/gc_content_genes.pkl"

with open(file_bacteria, "rb") as f:
    bacteria_gc_df = pickle.load(f)
    # 72690 rows - one column Bacteria_ID, one column GC_content
with open(file_genes, "rb") as f:
    genes_gc_df = pickle.load(f)
    # 6048 rows

# Go through each gene
gene_names = "/storage/enyaa/FINAL/gene_names.txt"
gene_names_df = pd.read_csv(gene_names, header=None, names=["Gene_name"])

euclidean_distances_all = []
gc_ratio_all = []

for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): # tqdm - progress bar

    gene_gc = genes_gc_df.loc[genes_gc_df['Gene_name'] == gene_name, 'GC_content'].values[0]

    # ONLY COMPUTE GC-RATIO BETWEEN GENES AND ITS MATCHES

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    eu_path = f"/storage/enyaa/FINAL/KMER/FOR_GENE_LENGTH/euclidean_split_genes_500bp/euclidean_df_{gene_name}.pkl"
    if not os.path.exists(eu_path): # skip genes that are shorter than 500 bp because those files don't exist
        continue
    euclidean_gene_df = pd.read_pickle(eu_path)

    # filter to only include matching bacteria
    filtered_euclidean_gene_df = euclidean_gene_df[euclidean_gene_df["Match_status"] == "Match"]

    if filtered_euclidean_gene_df.empty: # skip genes with no matches
        continue

    filtered_gc_df = bacteria_gc_df.merge(
        filtered_euclidean_gene_df[['Bacteria_ID']], on='Bacteria_ID', how='inner'
    )
    bacteria_gc = filtered_gc_df['GC_content'].to_numpy()

    # Compute GC ratio
    ratio = np.round(gene_gc / bacteria_gc, 4) # array of the ratio for one gene vs all genomes

    euclidean_distances_all.extend(filtered_euclidean_gene_df["Euclidean_distance"].values) 
    gc_ratio_all.extend(ratio) 


df_plot = pd.DataFrame({
    'Euclidean_distance': euclidean_distances_all,
    'GC_ratio': gc_ratio_all
})

# Scatterplot:
plt.figure(figsize=(8, 6))
plt.scatter(data=df_plot, x='Euclidean_distance', y='GC_ratio', alpha=1, s=10, color="darkorange")
plt.xlabel("Length-adjusted 5mer score", fontsize=16)
plt.xlim(0.015, 0.11)
plt.ylabel("GC-ratio", fontsize=16)
plt.tick_params(axis='both', labelsize=14)
#plt.title("GC-ratio vs length-adjusted 5mer score for matching genes and genomes")
plt.tight_layout()
plt.grid(True)
plt.savefig('/home/enyaa/gene_genome/scatterplot_GC_ratio_length_adjusted.png') 
plt.close()


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot ratio (vs length-adjusted) created in: {total_time} minutes!")
 """

# SCATTERPLOT DIFFERENCE for matches ----------------
# Calculates the GC-difference between genes and genomes - only for filtered matches
# Makes a scatterplot gc-difference vs length-adjusted 5mer score

start_time = time.time()

# Load the files with GC-content for genes and genomes
file_bacteria = "/storage/enyaa/FINAL/GC/gc_content_bacteria.pkl"
file_genes = "/storage/enyaa/FINAL/GC/gc_content_genes.pkl"

with open(file_bacteria, "rb") as f:
    bacteria_gc_df = pickle.load(f)
    # 72690 rows - one column Bacteria_ID, one column GC_content
with open(file_genes, "rb") as f:
    genes_gc_df = pickle.load(f)
    # 6048 rows

# Go through each gene
gene_names = "/storage/enyaa/FINAL/gene_names.txt"
gene_names_df = pd.read_csv(gene_names, header=None, names=["Gene_name"])

euclidean_distances_all = []
gc_diff_all = []

for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): # tqdm - progress bar

    gene_gc = genes_gc_df.loc[genes_gc_df['Gene_name'] == gene_name, 'GC_content'].values[0]

    # ONLY COMPUTE GC-DIFFERENCE BETWEEN GENES AND ITS MATCHES

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    eu_path = f"/storage/enyaa/FINAL/KMER/FOR_GENE_LENGTH/euclidean_split_genes_500bp/euclidean_df_{gene_name}.pkl"
    if not os.path.exists(eu_path): # skip genes that are shorter than 500 bp because those files don't exist
        continue
    euclidean_gene_df = pd.read_pickle(eu_path)

    # filter to only include matching bacteria
    filtered_euclidean_gene_df = euclidean_gene_df[euclidean_gene_df["Match_status"] == "Match"]

    if filtered_euclidean_gene_df.empty: # skip genes with no matches
        continue

    filtered_gc_df = bacteria_gc_df.merge(
        filtered_euclidean_gene_df[['Bacteria_ID']], on='Bacteria_ID', how='inner'
    )
    bacteria_gc = filtered_gc_df['GC_content'].to_numpy()

    # Compute GC difference
    diff = np.round(gene_gc - bacteria_gc, 4) # array of the difference for one gene vs all genomes

    euclidean_distances_all.extend(filtered_euclidean_gene_df["Euclidean_distance"].values) 
    gc_diff_all.extend(diff) 


df_plot = pd.DataFrame({
    'Euclidean_distance': euclidean_distances_all,
    'GC_diff': gc_diff_all
})

# Scatterplot:
plt.figure(figsize=(8, 6))
plt.scatter(data=df_plot, x='Euclidean_distance', y='GC_diff', alpha=1, s=10, color="darkorange")
plt.xlabel("Length-adjusted 5mer score", fontsize=16)
plt.xlim(0.015, 0.11)
plt.ylabel("GC-difference", fontsize=16)
plt.tick_params(axis='both', labelsize=14)
#plt.title("GC-diff vs length-adjusted 5mer score for matching genes and genomes")
plt.tight_layout()
plt.grid(True)
plt.savefig('/home/enyaa/gene_genome/scatterplot_GC_diff_length_adjusted.png') 
plt.close()


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot diff (vs length-adjusted) created in: {total_time} minutes!")

