
from Bio.SeqUtils import gc_fraction
import pandas as pd
import gzip
from Bio import SeqIO
import os
import pickle
import time
import matplotlib.pyplot as plt


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

# Load the files with GC-content for genes and genomes
file_bacteria = "/storage/enyaa/REVISED/GC/gc_content_bacteria.pkl"
file_genes = "/storage/enyaa/REVISED/GC/gc_content_genes.pkl"

with open(file_bacteria, "rb") as f:
    bacteria_gc_df = pickle.load(f)
with open(file_genes, "rb") as f:
    genes_gc_df = pickle.load(f)

# Ratio between gene GC content and bacteria GC content (gene/bacteria)
for gene in genes_gc_df.itertuples():
    gene_name = gene.Gene_name
    gene_gc = gene.GC_content

    # Create df for current gene (values = NaN)
    gene_df = pd.DataFrame(index=[gene_name], columns=bacteria_gc_df['Bacteria_ID'])
    gene_df.columns.name = None # take away the label "Bacteria_ID"

    for bacteria in bacteria_gc_df.itertuples():
        bacteria_id = bacteria.Bacteria_ID
        bacteria_gc = bacteria.GC_content

        # Calculate the ratio
        ratio = gene_gc / bacteria_gc

        # Store the ratio in the result DataFrame
        gene_df.at[gene_name, bacteria_id] = ratio

    # Save result to pickle - one for each gene
    output_dir = "/storage/enyaa/REVISED/GC/gc_ratio_split_genes"

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    gene_path = os.path.join(output_dir, f"{gene_name}.pkl")
    gene_df.to_pickle(gene_path)

print(f"All GC-files created!")
'''


# SCATTERPLOT for matches ------------------------------------
#'''
gene_names = "/home/enyaa/gene_genome/gene_names.txt"
gc_ratio_path = "/storage/enyaa/REVISED/GC/gc_ratio_split_genes/"

gene_names_df = pd.read_csv(gene_names, header=None, names=["Gene_name"])

euclidean_distances_all = []
gc_ratio_all = []

for gene_name in gene_names_df["Gene_name"][:3]:
    try:
        if "/" in gene_name:
            gene_name = gene_name.replace("/", "?")

        euclidean_gene_df = pd.read_pickle(f"/storage/enyaa/REVISED/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl")
        ratio_gene_df = pd.read_pickle(f"{gc_ratio_path}{gene_name}.pkl")

        # find matching bacteria
        path = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/taxonomy_results_{gene_name}.csv" # this contains matches
        if os.path.exists(path):
            taxonomy_gene_df = pd.read_csv(path, sep=",")
            matching_bacteria = taxonomy_gene_df[['Bacteria_ID']].drop_duplicates().tolist()
        else: # if taxonomy file doesn't exist - there are no matches for the gene, skip since we only are interested in matches in this code
            print(f"File not found: {path}")  
            continue  # Skip to the next gene
            
        # filter to only include matching bacteria
        filtered_euclidean_gene_df = euclidean_gene_df.loc[gene_name, matching_bacteria] # first column with Bacteria_ID and second column with Euclidean_distance - but no header!
        filtered_ratio_gene_df = ratio_gene_df.loc[gene_name, matching_bacteria] # first column with Bacteria_ID and second column with GC_ratio - but no header!

        euclidean_distances_all.extend(filtered_euclidean_gene_df.values) # list with euclidean distances
        gc_ratio_all.extend(filtered_ratio_gene_df.values) # list with gc_ratio, in same bacterial order as the euclidean_distances_all
    except Exception as e:
        print(f"Skipping {gene_name} due to error: {e}")

# Scatterplot:
plt.figure(figsize=(8, 6))
plt.scatter(euclidean_distances_all, gc_ratio_all, alpha=1, s=10)
plt.xlabel("Euclidean distance")
plt.ylabel("GC-ratio")
plt.title("GC-ratio vs euclidean distance for matching genes and genomes")
plt.grid(True)
plt.savefig('/home/enyaa/gene_genome/scatterplot_GC.png') 
plt.close()
#'''