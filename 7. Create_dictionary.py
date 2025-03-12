# Beräkna kmer distributions för alla genom genom att normalize the counts

import os
import pandas as pd
import pickle

'''
directory= "/storage/enyaa/REVISED/KMER/gene_kmer"
kmer_dict = {} # empty dictionary

for filename in os.listdir(directory): # goes through filenames in directory
    filepath = os.path.join(directory, filename) # creates the filepath
    gene_name = filename.split("_kmer_counts.txt")[0] # gets the gene name
    
    if "?" in gene_name:
        gene_name = gene_name.replace("?", "/") # ? = /
        #print(gene_name)
    df = pd.read_csv(filepath, sep="\t", header=None, names=["kmer", "count"])
    
    # Normalize counts:
    total = df["count"].sum()
    if total > 0: # Makes sure we do not divide by zero
        df["count"] /= total
    
    kmer_dict[gene_name] = df.set_index("kmer")["count"].to_dict() # Adds the normalized counts to dictionary


# Save dictionary to a pickle-file:
save_path = "/storage/enyaa/REVISED/KMER/gene_dist/gene_kmer_distributions.pkl"    
with open(save_path, "wb") as f:
    pickle.dump(kmer_dict, f)

'''

directory= "/storage/enyaa/REVISED/KMER/genome_kmer" 
kmer_dict = {}


for filename in os.listdir(directory):
    filepath = os.path.join(directory, filename)
    bacteria_id = filename.split("_kmc_counts.txt")[0]
    df = pd.read_csv(filepath, sep="\t", header=None, names=["kmer", "count"])
    
    # Normalize counts:
    total = df["count"].sum()
    if total > 0: # Makes sure we do not divide by zero
        df["count"] /= total
    
    kmer_dict[bacteria_id] = df.set_index("kmer")["count"].to_dict() # Adds the normalized counts to dictionary

# Save dictionary to a pickle-file:
save_path = "/storage/enyaa/REVISED/KMER/genome_dist/genome_kmer_distributions_11.pkl"    
with open(save_path, "wb") as f:
    pickle.dump(kmer_dict, f)
    
