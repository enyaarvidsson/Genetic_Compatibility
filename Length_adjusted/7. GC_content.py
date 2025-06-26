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
""" 
start_time = time.time()


# Load the files with GC-content for genes and genomes
file_bacteria = "/storage/enyaa/FINAL/GC/gc_content_bacteria.pkl"
file_genes = "/storage/enyaa/FINAL/GC/gc_content_genes.pkl"

with open(file_bacteria, "rb") as f:
    bacteria_gc_df = pickle.load(f)
    # should be 77182 rows - one column Bacteria_ID, one column GC_content
with open(file_genes, "rb") as f:
    genes_gc_df = pickle.load(f)

# Go through each gene
gene_names = "/storage/enyaa/FINAL/gene_names.txt"
gene_names_df = pd.read_csv(gene_names, header=None, names=["Gene_name"])

fivemer_score_all = []
gc_ratio_all = []
phylum_all = []

for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): # tqdm - progress bar

    gene_gc = genes_gc_df.loc[genes_gc_df['Gene_name'] == gene_name, 'GC_content'].values[0]

    # ONLY COMPUTE GC-RATIO BETWEEN GENES AND ITS MATCHES

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    fivemer_score_gene_df = pd.read_pickle(f"/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_{gene_name}.pkl")
    # borde finnas en kolumn med matchningar

    # filter to only include matching bacteria
    filtered_euclidean_gene_df = fivemer_score_gene_df[fivemer_score_gene_df["Match_status"] == "Match"]

    if filtered_euclidean_gene_df.empty: # skip genes with no matches
        continue

    filtered_gc_df = bacteria_gc_df.merge(
        filtered_euclidean_gene_df[['Bacteria_ID']], on='Bacteria_ID', how='inner'
    )
    bacteria_gc = filtered_gc_df['GC_content'].to_numpy()

    # Compute GC ratio
    ratio = np.round(gene_gc / bacteria_gc, 4) # array of the ratio for one gene vs all genomes

    phylum_all.extend(filtered_euclidean_gene_df["Phylum"].values)
    fivemer_score_all.extend(filtered_euclidean_gene_df["Euclidean_distance"].values) 
    gc_ratio_all.extend(ratio) 


df_plot = pd.DataFrame({
    'Euclidean_distance': fivemer_score_all,
    'GC_ratio': gc_ratio_all,
    'Phylum': phylum_all
})

# Scatterplot:
plt.figure(figsize=(8, 6))
plt.scatter(data=df_plot, x='Euclidean_distance', y='GC_ratio', alpha=1, s=10, color="darkorange")
plt.xlabel("5mer score", fontsize=16)
plt.xlim(0.015, 0.11)
plt.ylabel("GC-ratio", fontsize=16)
plt.tick_params(axis='both', labelsize=14)
#plt.title("GC-ratio vs euclidean distance for matching genes and genomes")
plt.tight_layout()
plt.grid(True)
plt.savefig('/home/enyaa/gene_genome/scatterplot_GC_ratio_euclidean.png') 
plt.close()


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot ratio created in: {total_time} minutes!")
 """

# SCATTERPLOT DIFFERENCE for matches ----------------

""" # KOPIERA OCH GÖR SAMMA SOM OVAN OM DET FUNKAR


# Calculates the GC-difference between genes and genomes - only for filtered matches
# Makes a scatterplot gc-diff vs euclidean
#'''

start_time = time.time()

# Load one euclidean distance file to take out the filtered bacteria
df = pd.read_pickle("/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_OXA-1095.pkl")
bacteria_ids = df['Bacteria_ID'].unique().tolist()
    # 77 182 bacteria_ids

# Load the files with GC-content for genes and genomes
file_bacteria = "/storage/enyaa/REVISED/GC/gc_content_bacteria.pkl"
file_genes = "/storage/enyaa/REVISED/GC/gc_content_genes.pkl"

with open(file_bacteria, "rb") as f:
    bacteria_gc_df = pickle.load(f)
with open(file_genes, "rb") as f:
    genes_gc_df = pickle.load(f)

bacteria_gc_filtered_df = bacteria_gc_df[bacteria_gc_df['Bacteria_ID'].isin(bacteria_ids)]
    # 77182 rows - one column Bacteria_ID, one column GC_content

# this makes code faster
available_taxonomy_files = set(os.listdir("/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/"))

# Go through each gene
gene_names = "/storage/enyaa/REVISED/gene_names.txt"
gene_names_df = pd.read_csv(gene_names, header=None, names=["Gene_name"])

euclidean_distances_500bp_all = []
gc_diff_all = []
phylum_all = []

for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): # tqdm - progress bar

    gene_gc = genes_gc_df.loc[genes_gc_df['Gene_name'] == gene_name, 'GC_content'].values[0]

    # ONLY COMPUTE GC-DIFF BETWEEN GENES AND ITS MATCHES

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    #euclidean_gene_df = pd.read_pickle(f"/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_{gene_name}.pkl")
    # for 500 bp:
    eu_path = f"/storage/enyaa/REVISED/KMER/FOR_GENE_LENGTH/euclidean_split_genes_500bp/euclidean_df_{gene_name}.pkl"
    if not os.path.exists(eu_path): # skip genes that are shorter than 500 bp because those files don't exist
        continue
    euclidean_gene_500bp_df = pd.read_pickle(eu_path).reset_index(drop=True).T.reset_index() # Switch to long format 
    euclidean_gene_500bp_df.columns = ['Bacteria_ID', 'Euclidean_distance']


    filename = f"taxonomy_results_{gene_name}.csv"
    if filename not in available_taxonomy_files:
        #print(f"File not found: {filename}")
        continue

    taxonomy_path = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/{filename}"
    taxonomy_gene_df = pd.read_csv(taxonomy_path)
    matching_df = taxonomy_gene_df[['Bacteria_ID', 'Phylum']]

    top_phyla = ['Pseudomonadota', 'Actinomycetota', 'Bacillota', 'Bacteroidota', 'Campylobacterota', 'Cyanobacteriota']
    matching_df = matching_df[matching_df['Phylum'].isin(top_phyla)]

    # filter to only include matching bacteria
    filtered_euclidean_gene_500bp_df = euclidean_gene_500bp_df.merge(matching_df, on='Bacteria_ID', how='inner')

    if filtered_euclidean_gene_500bp_df.empty: 
        #print("No matches for gene:", gene_name)
        continue

    filtered_gc_df = bacteria_gc_filtered_df.merge(
        filtered_euclidean_gene_500bp_df[['Bacteria_ID']], on='Bacteria_ID', how='inner'
    )
    bacteria_gc = filtered_gc_df['GC_content'].to_numpy()

    # Compute GC difference 
    diff = np.round(gene_gc - bacteria_gc, 4) # array of the difference for one gene vs all genomes

    phylum_all.extend(filtered_euclidean_gene_500bp_df["Phylum"].values)
    euclidean_distances_500bp_all.extend(filtered_euclidean_gene_500bp_df["Euclidean_distance"].values) 
    gc_diff_all.extend(diff) 


df_plot = pd.DataFrame({
    'Euclidean_distance': euclidean_distances_500bp_all,
    'GC_difference': gc_diff_all,
    'Phylum': phylum_all
})

# Scatterplot:
plt.figure(figsize=(8, 6))
#sns.scatterplot(data=df_plot, x='Euclidean_distance', y='GC_ratio', s=10)
plt.scatter(data=df_plot, x='Euclidean_distance', y='GC_difference', alpha=1, s=10, color="darkorange")
#plt.scatter(euclidean_distances_all, gc_diff_all, alpha=1, s=10)
plt.xlabel("Length-adjusted 5mer score", fontsize=16)
plt.xlim(0.015, 0.11)
plt.ylabel("GC-difference", fontsize=16)
plt.tick_params(axis='both', labelsize=14)
#plt.title("GC-ratio vs euclidean distance for matching genes and genomes")
#plt.legend(title='Phylum', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.grid(True)
plt.savefig('/home/enyaa/gene_genome/scatterplot_GC_diff_euclidean500bp.png') 
plt.close()


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot diff created in: {total_time} minutes!")
 """

