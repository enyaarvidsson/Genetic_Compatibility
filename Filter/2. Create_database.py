# Create BLAST database with our filtered bacteria

import pandas as pd
import os
import subprocess
from tqdm import tqdm
import gzip
from Bio import SeqIO


# Filtered bacteria ------------------
filtered_bacteria_df = pd.read_csv("/storage/enyaa/FINAL/filtered_bacteria.csv")
    # 72690 bacteria

# Get the file paths -----------------
file_paths_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_filepaths.tsv"
df = pd.read_csv(file_paths_file, sep="\t", header=None, names=["Filepath"]) 
df["Filename"] = df["Filepath"].apply(os.path.basename) # gets the filename (not whole filepath) 
df["Bacteria_ID"] = df["Filename"].str.split("_").str[:2].str.join("_") # gets the bacterial ID
df = df.drop(columns=["Filename"])
    # 1631873 bacteria

filtered_bacteria_df = pd.merge(filtered_bacteria_df[['Bacteria_ID']], df, on='Bacteria_ID', how='left')
    # 72690 bacteria, with filepaths
file_paths = filtered_bacteria_df['Filepath'].tolist()

# Create BLAST database --------------
# Open bacterial genome zip-files and create a fasta file
output_fasta = "/storage/enyaa/FINAL/BLAST/genomes.fasta"

with open(output_fasta, "w") as fasta_out:
    
    for path in tqdm(file_paths, desc="Processing bacteria"): # goes through paths to the bacterial genome zip-files
        try:
            filename = os.path.basename(path) # gets the filename (not whole filepath) 
            bacteria_id = "_".join(filename.split("_")[:2]) # gets the bacterial ID
            
            with gzip.open(path, "rt", encoding="utf-8") as gz_file: # opens the zip file for reading
                for record in SeqIO.parse(gz_file, "fasta"):
        
                    record.id = bacteria_id + "|" + record.id  # adds bacteria_id to header in the fasta file
                    
                    SeqIO.write(record, fasta_out, "fasta") # writes in the output file
                    
        except Exception as e:
            print(f"Error opening {path}: {e}")      


# Create BLAST database:
blast_path = "/storage/enyaa/FINAL/BLAST/DB/blast" # name for BLAST database
subprocess.run(["makeblastdb", "-in", output_fasta, "-dbtype", "nucl", "-out", blast_path]) # runs this command in the terminal to create a database


print(f"Database created!")
