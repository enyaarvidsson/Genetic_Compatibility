# This script computes k-mer (5-mer) distributions
    # for the genes
    # and for the genomes
# by normalizing the counts.

import os
import pandas as pd
import pickle
from tqdm import tqdm


# FOR THE GENES --------------------

""" directory= "/storage/enyaa/FINAL/KMER/gene_kmer"

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
save_path = "/storage/enyaa/FINAL/KMER/gene_kmer_distributions.pkl"     
with open(save_path, "wb") as f:
    pickle.dump(kmer_dict, f) 
    
print("Created a pkl-file for the genes kmer distributions!")
"""


# FOR THE GENOMES -----------------

directory= "/storage/enyaa/FINAL/KMER/genome_kmer" 
kmer_dict = {}

for filename in tqdm(os.listdir(directory), desc="Processing genomes"):
    filepath = os.path.join(directory, filename)
    bacteria_id = filename.split("_kmc_counts.txt")[0]
    df = pd.read_csv(filepath, sep="\t", header=None, names=["kmer", "count"])
    
    # Normalize counts:
    total = df["count"].sum()
    if total > 0: # Makes sure we do not divide by zero
        df["count"] /= total
    
    kmer_dict[bacteria_id] = df.set_index("kmer")["count"].to_dict() # Adds the normalized counts to dictionary

# Save dictionary to a pickle-file:
save_path = "/storage/enyaa/FINAL/KMER/genome_kmer_distributions.pkl"    
with open(save_path, "wb") as f:
    pickle.dump(kmer_dict, f)
    
print("Created a pkl-file for the genomes kmer distributions!")