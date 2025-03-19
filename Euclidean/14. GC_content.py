
from Bio.SeqUtils import gc_fraction
import pandas as pd
import gzip
from Bio import SeqIO
import os


# GC-CONTENT FOR BACTERIA ------------------------------------

file_paths_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_filepaths.tsv"
df = pd.read_csv(file_paths_file, sep="\t", header=None) 
file_paths = df[0].tolist()

results = []

for path in file_paths:

    filename = os.path.basename(path) # gets the filename (not whole filepath) 
    bacteria_id = "_".join(filename.split("_")[:2]) # gets the bacterial ID

    # all sequences in one fasta file:
    all_sequences = ""  # empty string

    with gzip.open(path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            all_sequences += str(record.seq)  # merge all sequences in one fasta-file

    # Compute GC-content for the entire dataset
    gc_content = round(gc_fraction(all_sequences) * 100, 2)  # Convert fraction to percentage

    #print(f"Total GC-content for {bacteria_id}: {gc_content:.2f}%")

    results.append([bacteria_id, gc_content])

# Save to csv-file
output_file = "/storage/enyaa/REVISED/GC/gc_content_bacteria.csv"
df_output = pd.DataFrame(results, columns=["Bacteria_ID", "GC_content"])
df_output.to_csv(output_file, index=False)



# GC-CONTENT FOR GENES ---------------------------------------
'''
file_path = "/storage/enyaa/nucleotide_fasta_protein_homolog_model.fasta"
output_file = "gc_content_genes.csv"

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
df_output.to_csv(output_file, index=False, float_format="%.2f")
'''



# MAKE DF FOR RATIO between genes and genomes ----------------

# !!!!!!!!!!!!!!!!!!!!!!!!
# DEN HÄR KODEN ÄR BARA INLAGD FRÅN CHATGPT, INTE TESTAD
# !!!!!!!!!!!!!!!!!!!!!!!!

import pandas as pd

# Step 1: Load the CSV files into pandas DataFrames
genes_df = pd.read_csv('gc_content_gene.csv')
bacteria_df = pd.read_csv('gc_content_bacteria.csv')

# Step 2: Create a DataFrame to store the results (ratio of gene GC content to bacteria GC content)
# Initialize an empty DataFrame with genes as rows and bacteria as columns
result_df = pd.DataFrame(index=genes_df['Gene_name'], columns=bacteria_df['Bacteria_ID'])

# Step 3: Calculate the ratio between gene GC content and bacteria GC content
for gene in genes_df.itertuples():
    for bacteria in bacteria_df.itertuples():
        gene_name = gene.Gene_name
        bacteria_id = bacteria.Bacteria_ID
        gene_gc = gene.GC_content
        bacteria_gc = bacteria.GC_content
        # Calculate the ratio
        ratio = gene_gc / bacteria_gc
        # Store the ratio in the result DataFrame
        result_df.at[gene_name, bacteria_id] = ratio

# Step 4: Save the result DataFrame to a CSV file
result_df.to_csv('gc_content_ratios.csv')

# Optional: Display the result
print(result_df)



