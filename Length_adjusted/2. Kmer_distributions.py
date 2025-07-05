# This script computes k-mer (5-mer) distributions
# for the length-adjusted genes, by normalizing the counts.

import os
import pandas as pd
import pickle
from tqdm import tqdm


directory= "./FINAL/KMER/FOR_GENE_LENGTH/gene_kmer_500bp"
kmer_dict = {} # empty dictionary

for filename in tqdm(os.listdir(directory), desc="Processing genes"): # goes through filenames in directory
    filepath = os.path.join(directory, filename) # creates the filepath
    gene_name = filename.split("_kmer_counts.txt")[0] # gets the gene name
    
    if "?" in gene_name:
        gene_name = gene_name.replace("?", "/") # ? = /
    df = pd.read_csv(filepath, sep="\t", header=None, names=["kmer", "count"])
    
    # Normalize counts:
    total = df["count"].sum()
    if total > 0: # Makes sure we do not divide by zero
        df["count"] /= total
    
    kmer_dict[gene_name] = df.set_index("kmer")["count"].to_dict() # Adds the normalized counts to dictionary


# Save dictionary to a pickle-file:  
save_path = "./FINAL/KMER/FOR_GENE_LENGTH/gene_kmer_distributions_500bp.pkl"  
with open(save_path, "wb") as f:
    pickle.dump(kmer_dict, f)

print("Created the kmer-distribution pkl-file, for the length-adjusted genes!")
