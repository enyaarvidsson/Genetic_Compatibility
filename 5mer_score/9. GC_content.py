# This script includes four parts:
    # Calculating GC-content for the bacteria, and storing it in a pkl file
    # Calculating GC-content for the genes, and storing it in a pkl file
    # Calculating the GC-ratio and making a scatterplot (GC-ratio vs 5mer score)
    # Calculating the GC-difference and making a scatterplot (GC-difference vs 5mer score)

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


# GC-CONTENT FOR BACTERIA ------------------------------------
# GC-content for the bacteria

""" start_time = time.time()

# Filtered bacteria 
filtered_bacteria_df = pd.read_csv("./FINAL/filtered_bacteria.csv")
    # 72690 bacteria

# Get the file paths 
file_paths_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_filepaths.tsv"
df = pd.read_csv(file_paths_file, sep="\t", header=None, names=["Filepath"]) 
df["Filename"] = df["Filepath"].apply(os.path.basename) # gets the filename (not whole filepath) 
df["Bacteria_ID"] = df["Filename"].str.split("_").str[:2].str.join("_") # gets the bacterial ID
df = df.drop(columns=["Filename"])
    # 1631873 bacteria

filtered_bacteria_df = pd.merge(filtered_bacteria_df[['Bacteria_ID']], df, on='Bacteria_ID', how='left')
    # 72690 bacteria, with filepaths
file_paths = filtered_bacteria_df['Filepath'].tolist()

results = []
for path in tqdm(file_paths, desc="Processing bacteria"):

    filename = os.path.basename(path) # gets the filename (not whole filepath) 
    bacteria_id = "_".join(filename.split("_")[:2]) # gets the bacterial ID

    sequences = [] 

    with gzip.open(path, "rt") as handle:
        sequences = [str(record.seq) for record in SeqIO.parse(handle, "fasta")] # store all sequences in a fasta file, in a list
    
    all_sequences = "".join(sequences) # join all sequences    
        
    # Compute GC-content for the entire dataset
    gc_content = round(gc_fraction(all_sequences) * 100, 2)  # Convert fraction to percentage

    results.append([bacteria_id, gc_content])

# Save to pickle
df_output = pd.DataFrame(results, columns=["Bacteria_ID", "GC_content"])
directory = os.path.join('.', 'FINAL', 'GC')
os.makedirs(directory, exist_ok=True) 
output_file = "./FINAL/GC/gc_content_bacteria.pkl"
with open(output_file, "wb") as f:
    pickle.dump(df_output, f)


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Bacteria gc-content file created in: {total_time} minutes") """


# GC-CONTENT FOR GENES ---------------------------------------
""" file_path = "./nucleotide_fasta_protein_homolog_model.fasta"
output_file = "./FINAL/GC/gc_content_genes.pkl"

results = []

for record in tqdm(SeqIO.parse(file_path, "fasta"), desc="Processing genes"):
    gene = record.id.split(" ")[0] 
    gene_name = gene.split("|")[5]

    # compute GC-content
    gc_content = round(gc_fraction(record.seq) * 100, 2)

    results.append([gene_name, gc_content])

# Save to pkl
df_output = pd.DataFrame(results, columns=["Gene_name", "GC_content"])

with open(output_file, "wb") as f:
    pickle.dump(df_output, f)

print(f"Saved pickle results to {output_file}")
 """

# SCATTERPLOT RATIO for matches ------------------------------
# Calculates the GC-ratio between genes and genomes - only for filtered matches
# Makes a scatterplot gc-ratio vs 5mer score
""" 
start_time = time.time()

# Load the files with GC-content for genes and genomes
file_bacteria = "./FINAL/GC/gc_content_bacteria.pkl"
file_genes = "./FINAL/GC/gc_content_genes.pkl"

with open(file_bacteria, "rb") as f:
    bacteria_gc_df = pickle.load(f)
    # 72690 rows - one column Bacteria_ID, one column GC_content
with open(file_genes, "rb") as f:
    genes_gc_df = pickle.load(f)
    # 6048 rows

# Go through each gene
gene_names = "./FINAL/gene_names.txt"
gene_names_df = pd.read_csv(gene_names, header=None, names=["Gene_name"])

euclidean_distances_all = []
gc_ratio_all = []

for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): # tqdm - progress bar

    gene_gc = genes_gc_df.loc[genes_gc_df['Gene_name'] == gene_name, 'GC_content'].values[0]

    # ONLY COMPUTE GC-RATIO BETWEEN GENES AND ITS MATCHES

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    euclidean_gene_df = pd.read_pickle(f"./FINAL/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl")

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
plt.xlabel("5mer score", fontsize=16)
plt.xlim(0.015, 0.11)
plt.ylabel("GC-ratio", fontsize=16)
plt.tick_params(axis='both', labelsize=14)
#plt.title("GC-ratio vs euclidean distance for matching genes and genomes")
plt.tight_layout()
plt.grid(True)
plt.savefig('./FINAL/GC/scatterplot_GC_ratio_5mer_score.png') 
plt.close()


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot ratio (vs 5mer score) created in: {total_time} minutes!")
 """

# SCATTERPLOT DIFFERENCE for matches ----------------
# Calculates the GC-difference between genes and genomes - only for filtered matches
# Makes a scatterplot gc-difference vs 5mer score

start_time = time.time()

# Load the files with GC-content for genes and genomes
file_bacteria = "./FINAL/GC/gc_content_bacteria.pkl"
file_genes = "./FINAL/GC/gc_content_genes.pkl"

with open(file_bacteria, "rb") as f:
    bacteria_gc_df = pickle.load(f)
    # 72690 rows - one column Bacteria_ID, one column GC_content
with open(file_genes, "rb") as f:
    genes_gc_df = pickle.load(f)
    # 6048 rows

# Go through each gene
gene_names = "./FINAL/gene_names.txt"
gene_names_df = pd.read_csv(gene_names, header=None, names=["Gene_name"])

euclidean_distances_all = []
gc_diff_all = []

for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): # tqdm - progress bar

    gene_gc = genes_gc_df.loc[genes_gc_df['Gene_name'] == gene_name, 'GC_content'].values[0]

    # ONLY COMPUTE GC-DIFFERENCE BETWEEN GENES AND ITS MATCHES

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    euclidean_gene_df = pd.read_pickle(f"./FINAL/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl")

    # filter to only include matching bacteria
    filtered_euclidean_gene_df = euclidean_gene_df[euclidean_gene_df["Match_status"] == "Match"]

    if filtered_euclidean_gene_df.empty: # skip genes with no matches
        continue

    filtered_gc_df = bacteria_gc_df.merge(
        filtered_euclidean_gene_df[['Bacteria_ID']], on='Bacteria_ID', how='inner'
    )
    bacteria_gc = filtered_gc_df['GC_content'].to_numpy()

    # Compute GC difference
    diff = np.round(gene_gc - bacteria_gc, 4) # array of the diff for one gene vs all genomes

    euclidean_distances_all.extend(filtered_euclidean_gene_df["Euclidean_distance"].values) 
    gc_diff_all.extend(diff) 


df_plot = pd.DataFrame({
    'Euclidean_distance': euclidean_distances_all,
    'GC_diff': gc_diff_all
})

# Scatterplot:
plt.figure(figsize=(8, 6))
plt.scatter(data=df_plot, x='Euclidean_distance', y='GC_diff', alpha=1, s=10, color="darkorange")
plt.xlabel("5mer score", fontsize=16)
plt.xlim(0.015, 0.11)
plt.ylabel("GC-difference", fontsize=16)
plt.tick_params(axis='both', labelsize=14)
#plt.title("GC-difference vs euclidean distance for matching genes and genomes")
plt.tight_layout()
plt.grid(True)
plt.savefig('./FINAL/GC/scatterplot_GC_diff_5mer_score.png') 
plt.close()


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot difference (vs 5mer score) created in: {total_time} minutes!")

