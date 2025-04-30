
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


tRNA_score = "tRNA_score_one_sided"


start_time = time.time()

# Load one tRNA_score file to take out the filtered bacteria 
filepath_tRNA = "/storage/jolunds/REVISED/tRNA/tRNA_score_new/tRNA_score_tet(Q).csv"
df = pd.read_csv(filepath_tRNA) 

# only top phyla
top_phyla = df["Phylum"].value_counts().head(6)
df = df[df["Phylum"].isin(top_phyla.index)]
bacteria_ids = df['Bacteria_ID'].unique().tolist()
    # 67 089 bacteria_ids

# Load the files with GC-content for genes and genomes 
file_bacteria = "/storage/enyaa/REVISED/GC/gc_content_bacteria.pkl"
file_genes = "/storage/enyaa/REVISED/GC/gc_content_genes.pkl"

with open(file_bacteria, "rb") as f:
    bacteria_gc_df = pickle.load(f)
with open(file_genes, "rb") as f:
    genes_gc_df = pickle.load(f)

bacteria_gc_filtered_df = bacteria_gc_df[bacteria_gc_df['Bacteria_ID'].isin(bacteria_ids)]
    # 67089 rows - one column Bacteria_ID, one column GC_content

# this makes code faster
available_taxonomy_files = set(os.listdir("/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/"))

# Go through each gene
gene_names = "/storage/enyaa/REVISED/gene_names.txt"
gene_names_df = pd.read_csv(gene_names, header=None, names=["Gene_name"])

tRNA_score_all = []
gc_ratio_all = []

for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): 

    gene_gc = genes_gc_df.loc[genes_gc_df['Gene_name'] == gene_name, 'GC_content'].values[0]

    # ONLY COMPUTE GC-RATIO BETWEEN GENES AND ITS MATCHES

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    tRNA_gene_df = pd.read_csv(f"/storage/jolunds/REVISED/tRNA/tRNA_score_new/tRNA_score_{gene_name}.csv") 
        # 70911 rows 

    filename = f"taxonomy_results_{gene_name}.csv"
    if filename not in available_taxonomy_files:
        #print(f"File not found: {filename}")
        continue

    taxonomy_path = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/{filename}"
    taxonomy_gene_df = pd.read_csv(taxonomy_path)
    matching_df = taxonomy_gene_df[['Bacteria_ID', 'Phylum']]

    # filter to only include matching bacteria
    filtered_tRNA_gene_df = tRNA_gene_df.merge(matching_df, on='Bacteria_ID', how='inner')

    if filtered_tRNA_gene_df.empty: 
        #print("No matches for gene:", gene_name)
        continue

    filtered_gc_df = bacteria_gc_filtered_df.merge(
        filtered_tRNA_gene_df[['Bacteria_ID']], on='Bacteria_ID', how='inner'
    )
    if filtered_gc_df.empty:
        #print("här är felet")
        continue
    bacteria_gc = filtered_gc_df['GC_content'].to_numpy()

    # Compute GC ratio
    ratio = np.round(gene_gc / bacteria_gc, 4) # array of the ratio for one gene vs all genomes

    tRNA_score_all.extend(filtered_tRNA_gene_df[tRNA_score].values) 
    gc_ratio_all.extend(ratio)

df_plot = pd.DataFrame({
    'tRNA_score': tRNA_score_all,
    'GC_ratio': gc_ratio_all
    #'Phylum': phylum_all
})

# Scatterplot:
plt.figure(figsize=(8, 6))
plt.scatter(data=df_plot, x='tRNA_score', y='GC_ratio', alpha=1, s=10)
plt.xlabel("tRNA score", fontsize=16)
plt.ylabel("GC-ratio", fontsize=16)
plt.tick_params(axis='both', labelsize=14)

if tRNA_score == "tRNA_score_one_sided":
    tRNA_score_title = "one-sided"
    tRNA_score_nr = "1"
else:
    tRNA_score_title = "two-sided"
    tRNA_score_nr = "2"

plt.title(f"GC-ratio vs tRNA score ({tRNA_score_title}) for matching genes and genomes")
#plt.legend(title='Phylum', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.grid(True)
plt.savefig(f'/home/enyaa/gene_genome/GC_ratio_vs_tRNA{tRNA_score_nr}.png') 
plt.close()


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot ratio created in: {total_time} minutes!")



# SCATTERPLOT DIFFERENCE for matches -------------------------
# Calculates the GC-difference between genes and genomes - only for filtered matches
# Makes a scatterplot gc-diff vs tRNA score

'''
tRNA_score = "tRNA_score_one_sided"


start_time = time.time()

# Load one tRNA_score file to take out the filtered bacteria 
filepath_tRNA = "/storage/jolunds/REVISED/tRNA/tRNA_score_new/tRNA_score_tet(Q).csv"
df = pd.read_csv(filepath_tRNA) 

# only top phyla
top_phyla = df["Phylum"].value_counts().head(6)
df = df[df["Phylum"].isin(top_phyla.index)]

bacteria_ids = df['Bacteria_ID'].unique().tolist()
    # 67 089 bacteria_ids

# Load the files with GC-content for genes and genomes 
file_bacteria = "/storage/enyaa/REVISED/GC/gc_content_bacteria.pkl"
file_genes = "/storage/enyaa/REVISED/GC/gc_content_genes.pkl"

with open(file_bacteria, "rb") as f:
    bacteria_gc_df = pickle.load(f)
with open(file_genes, "rb") as f:
    genes_gc_df = pickle.load(f)

bacteria_gc_filtered_df = bacteria_gc_df[bacteria_gc_df['Bacteria_ID'].isin(bacteria_ids)]
    # 67089 rows - one column Bacteria_ID, one column GC_content

# this makes code faster
available_taxonomy_files = set(os.listdir("/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/"))

# Go through each gene
gene_names = "/storage/enyaa/REVISED/gene_names.txt"
gene_names_df = pd.read_csv(gene_names, header=None, names=["Gene_name"])

tRNA_score_all = []
gc_diff_all = []

for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): 

    gene_gc = genes_gc_df.loc[genes_gc_df['Gene_name'] == gene_name, 'GC_content'].values[0]

    # ONLY COMPUTE GC-DIFF BETWEEN GENES AND ITS MATCHES

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    tRNA_score_df = pd.read_csv(f"/storage/jolunds/REVISED/tRNA/tRNA_score_new/tRNA_score_{gene_name}.csv") 
    
    filename = f"taxonomy_results_{gene_name}.csv"
    if filename not in available_taxonomy_files:
        #print(f"File not found: {filename}")
        continue

    taxonomy_path = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/{filename}"
    taxonomy_gene_df = pd.read_csv(taxonomy_path)
    matching_df = taxonomy_gene_df[['Bacteria_ID', 'Phylum']]

    # filter to only include matching bacteria
    filtered_tRNA_score_df = tRNA_score_df.merge(matching_df, on='Bacteria_ID', how='inner')

    if filtered_tRNA_score_df.empty: 
        #print("No matches for gene:", gene_name)
        continue

    filtered_gc_df = bacteria_gc_filtered_df.merge(
        filtered_tRNA_score_df[['Bacteria_ID']], on='Bacteria_ID', how='inner'
    )
    if filtered_gc_df.empty:
        #print("här är felet")
        continue
    bacteria_gc = filtered_gc_df['GC_content'].to_numpy()

    # Compute GC difference
    diff = np.round(gene_gc - bacteria_gc, 4) # array of the difference for one gene vs all genomes

    tRNA_score_all.extend(filtered_tRNA_score_df[tRNA_score].values) 
    gc_diff_all.extend(diff) 


df_plot = pd.DataFrame({
    'tRNA_score': tRNA_score_all,
    'GC_difference': gc_diff_all
    #'Phylum': phylum_all
})

# Scatterplot:
plt.figure(figsize=(8, 6))
plt.scatter(data=df_plot, x='tRNA_score', y='GC_difference', alpha=1, s=10)
plt.xlabel("tRNA score", fontsize=16)
plt.ylabel("GC-difference", fontsize=16)
plt.tick_params(axis='both', labelsize=14)

if tRNA_score == "tRNA_score_one_sided":
    tRNA_score_title = "one-sided"
    tRNA_score_nr = "1"
else:
    tRNA_score_title = "two-sided"
    tRNA_score_nr = "2"

plt.title(f"GC-difference vs tRNA score ({tRNA_score_title}) for matching genes and genomes")
#plt.legend(title='Phylum', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.grid(True)
plt.savefig(f'/home/enyaa/gene_genome/GC_diff_vs_tRNA{tRNA_score_nr}.png') 
plt.close()


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot difference created in: {total_time} minutes!")
'''

