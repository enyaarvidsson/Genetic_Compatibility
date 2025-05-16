
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
'''

start_time = time.time()

file_paths_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_filepaths.tsv"
df = pd.read_csv(file_paths_file, sep="\t", header=None) 
file_paths = df[0].tolist()

results = []

for path in file_paths:

    filename = os.path.basename(path) # gets the filename (not whole filepath) 
    bacteria_id = "_".join(filename.split("_")[:2]) # gets the bacterial ID

    sequences = [] 

    with gzip.open(path, "rt") as handle:
        sequences = [str(record.seq) for record in SeqIO.parse(handle, "fasta")] # store all sequences in a fasta file, in a list
    
    all_sequences = "".join(sequences) # join all sequences    
        
    # Compute GC-content for the entire dataset
    gc_content = round(gc_fraction(all_sequences) * 100, 2)  # Convert fraction to percentage

    #print(f"Total GC-content for {bacteria_id}: {gc_content:.2f}%")

    results.append([bacteria_id, gc_content])

# Save to pickle
df_output = pd.DataFrame(results, columns=["Bacteria_ID", "GC_content"])
output_file = "/storage/enyaa/REVISED/GC/gc_content_bacteria.pkl"
with open(output_file, "wb") as f:
    pickle.dump(df_output, f)


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Bacteria gc-content file created in: {total_time} minutes")
'''


# GC-CONTENT FOR GENES ---------------------------------------
'''
file_path = "/storage/enyaa/nucleotide_fasta_protein_homolog_model.fasta"
output_file = "/storage/enyaa/REVISED/GC/gc_content_genes.pkl"

results = []

for record in SeqIO.parse(file_path, "fasta"):
    gene = record.id.split(" ")[0] 
    gene_name = gene.split("|")[5]

    # compute GC-content
    gc_content = round(gc_fraction(record.seq) * 100, 2)
    
    #print(f"Gene: {gene_name}, GC-content: {gc_content:.2f}%")

    results.append([gene_name, gc_content])

# Save to CSV
df_output = pd.DataFrame(results, columns=["Gene_name", "GC_content"])
#df_output.to_csv(output_file, index=False, float_format="%.2f")
with open(output_file, "wb") as f:
    pickle.dump(df_output, f)

print(f"Saved pickle results to {output_file}")
'''

# MAKE DF FOR RATIO between genes and genomes ----------------
'''

start_time = time.time()

# Load the files with GC-content for genes and genomes
file_bacteria = "/storage/enyaa/REVISED/GC/gc_content_bacteria.pkl"
file_genes = "/storage/enyaa/REVISED/GC/gc_content_genes.pkl"

with open(file_bacteria, "rb") as f:
    bacteria_gc_df = pickle.load(f)
with open(file_genes, "rb") as f:
    genes_gc_df = pickle.load(f)

# Bacteria ids and GC content
bacteria_ids = bacteria_gc_df["Bacteria_ID"].values # array with all bacteria_ids
bacteria_gc = bacteria_gc_df["GC_content"].astype(np.float32).values # array with gc_content for bacteria

output_dir = "/storage/enyaa/REVISED/GC/gc_ratio_split_genes"

# Go through each gene
for _, row in genes_gc_df.iterrows(): # _ is the nr of the row starting at 0, row is the gene_name and gc_content for one gene
    gene_name = row["Gene_name"]
    gene_gc = row["GC_content"]

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    gene_path = os.path.join(output_dir, f"{gene_name}.csv")

    # Compute GC ratio 
    ratio = np.round(gene_gc / bacteria_gc, 4) # array of the ratio for one gene vs all genomes

    # Save one file for each gene 
    pd.DataFrame([ratio], columns=bacteria_ids).to_csv(gene_path, index=False, float_format="%.2f")


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"All GC-files created in: {total_time} minutes!")
'''


# SCATTERPLOT RATIO for matches ------------------------------------
'''

start_time = time.time()

gene_names = "/storage/enyaa/REVISED/gene_names.txt"
gc_ratio_path = "/storage/enyaa/REVISED/GC/gc_ratio_split_genes/"

gene_names_df = pd.read_csv(gene_names, header=None, names=["Gene_name"])

euclidean_distances_all = []
gc_ratio_all = []

for gene_name in gene_names_df["Gene_name"]:
    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    euclidean_gene_df = pd.read_pickle(f"/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_{gene_name}.pkl")
    ratio_gene_df = pd.read_csv(f"{gc_ratio_path}{gene_name}.csv")

    # find matching bacteria
    path = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/taxonomy_results_{gene_name}.csv" # this contains matches
    if os.path.exists(path):
        taxonomy_gene_df = pd.read_csv(path, sep=",")
        matching_bacteria = taxonomy_gene_df['Bacteria_ID'].drop_duplicates().tolist() # contains all the matching bacteria_id to the current gene
    else: # if taxonomy file doesn't exist - there are no matches for the gene, skip since we only are interested in matches in this code
        print(f"File not found: {path}")  
        continue  # Skip to the next gene
            
    # filter to only include matching bacteria
    filtered_euclidean_gene_df = euclidean_gene_df[euclidean_gene_df['Bacteria_ID'].isin(matching_bacteria)] # first column with Bacteria_ID and second column with Euclidean_distance
    
    if filtered_euclidean_gene_df.empty: 
        print("No matches for gene:", gene_name)
        continue
    
    # only matching bacteria, filtered
    matching_bacteria_filtered = filtered_euclidean_gene_df['Bacteria_ID'].tolist()

    filtered_ratio_gene = ratio_gene_df.loc[0, ratio_gene_df.columns.isin(matching_bacteria_filtered)] 
    
    # make into a df
    filtered_ratio_gene_df = pd.DataFrame({
        "Bacteria_ID": filtered_ratio_gene.index,
        "GC_ratio": filtered_ratio_gene.values
    })

    # merge euclidean and ratio
    merged_df = pd.merge(filtered_euclidean_gene_df, filtered_ratio_gene_df, on="Bacteria_ID", how="inner")

    euclidean_distances_all.extend(merged_df["Euclidean_distance"].values) 
    gc_ratio_all.extend(merged_df["GC_ratio"].values) 

    
# Scatterplot:
plt.figure(figsize=(8, 6))
plt.scatter(euclidean_distances_all, gc_ratio_all, alpha=1, s=10)
plt.xlabel("Euclidean distance")
plt.ylabel("GC-ratio")
plt.title("GC-ratio vs euclidean distance for matching genes and genomes")
plt.grid(True)
plt.savefig('/home/enyaa/gene_genome/scatterplot_GC.png') 
plt.close()

end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot created in: {total_time} minutes!")
'''


# SCATTERPLOT RATIO for matches - FASTER ----------------
# Calculates the GC-ratio between genes and genomes - only for filtered matches
# Makes a scatterplot gc-ratio vs euclidean
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
gc_ratio_all = []
phylum_all = []

for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): # tqdm - progress bar

    gene_gc = genes_gc_df.loc[genes_gc_df['Gene_name'] == gene_name, 'GC_content'].values[0]

    # ONLY COMPUTE GC-RATIO BETWEEN GENES AND ITS MATCHES

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    #euclidean_gene_500bp_df = pd.read_pickle(f"/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_{gene_name}.pkl")
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

    # Compute GC ratio
    ratio = np.round(gene_gc / bacteria_gc, 4) # array of the ratio for one gene vs all genomes

    phylum_all.extend(filtered_euclidean_gene_500bp_df["Phylum"].values)
    euclidean_distances_500bp_all.extend(filtered_euclidean_gene_500bp_df["Euclidean_distance"].values) 
    gc_ratio_all.extend(ratio) 


df_plot = pd.DataFrame({
    'Euclidean_distance': euclidean_distances_500bp_all,
    'GC_ratio': gc_ratio_all,
    'Phylum': phylum_all
})

# Scatterplot:
plt.figure(figsize=(8, 6))
#sns.scatterplot(data=df_plot, x='Euclidean_distance', y='GC_ratio', s=10)
plt.scatter(data=df_plot, x='Euclidean_distance', y='GC_ratio', alpha=1, s=10, color="darkorange")
#plt.scatter(euclidean_distances_all, gc_diff_all, alpha=1, s=10)
plt.xlabel("Length-adjusted 5mer score", fontsize=16)
plt.xlim(0.015, 0.11)
plt.ylabel("GC-ratio", fontsize=16)
plt.tick_params(axis='both', labelsize=14)
#plt.title("GC-ratio vs euclidean distance for matching genes and genomes")
#plt.legend(title='Phylum', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.grid(True)
plt.savefig('/home/enyaa/gene_genome/scatterplot_GC_ratio_euclidean_500bp.png') 
plt.close()


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot ratio created in: {total_time} minutes!")
#'''





# SCATTERPLOT DIFFERENCE for matches ----------------
# Calculates the GC-difference between genes and genomes - only for filtered matches
# Makes a scatterplot gc-diff vs euclidean
'''

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

euclidean_distances_all = []
gc_diff_all = []
phylum_all = []

for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): # tqdm - progress bar

    gene_gc = genes_gc_df.loc[genes_gc_df['Gene_name'] == gene_name, 'GC_content'].values[0]

    # ONLY COMPUTE GC-DIFF BETWEEN GENES AND ITS MATCHES

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    euclidean_gene_df = pd.read_pickle(f"/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_{gene_name}.pkl")

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
    filtered_euclidean_gene_df = euclidean_gene_df.merge(matching_df, on='Bacteria_ID', how='inner')

    if filtered_euclidean_gene_df.empty: 
        #print("No matches for gene:", gene_name)
        continue

    filtered_gc_df = bacteria_gc_filtered_df.merge(
        filtered_euclidean_gene_df[['Bacteria_ID']], on='Bacteria_ID', how='inner'
    )
    bacteria_gc = filtered_gc_df['GC_content'].to_numpy()

    # Compute GC difference 
    diff = np.round(gene_gc - bacteria_gc, 4) # array of the difference for one gene vs all genomes

    phylum_all.extend(filtered_euclidean_gene_df["Phylum"].values)
    euclidean_distances_all.extend(filtered_euclidean_gene_df["Euclidean_distance"].values) 
    gc_diff_all.extend(diff) 


df_plot = pd.DataFrame({
    'Euclidean_distance': euclidean_distances_all,
    'GC_difference': gc_diff_all,
    'Phylum': phylum_all
})

# Scatterplot:
plt.figure(figsize=(8, 6))
#sns.scatterplot(data=df_plot, x='Euclidean_distance', y='GC_ratio', s=10)
plt.scatter(data=df_plot, x='Euclidean_distance', y='GC_difference', alpha=1, s=10, color="darkorange")
#plt.scatter(euclidean_distances_all, gc_diff_all, alpha=1, s=10)
plt.xlabel("5mer score", fontsize=16)
plt.ylabel("GC-difference", fontsize=16)
plt.tick_params(axis='both', labelsize=14)
#plt.title("GC-ratio vs euclidean distance for matching genes and genomes")
#plt.legend(title='Phylum', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.grid(True)
plt.savefig('/home/enyaa/gene_genome/scatterplot_GC_diff_euclidean.png') 
plt.close()


end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Scatterplot diff created in: {total_time} minutes!")
'''

