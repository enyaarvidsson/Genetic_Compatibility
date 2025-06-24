# FLYTTAD!

# Create database

import gzip
from Bio import SeqIO
import csv
import subprocess
import os
import time
import pandas as pd

start_time = time.time()

'''
# Bacterial genomes:
# get the first 100k lines in the file - in a more compicated way than below
filepaths = "/storage/shared/data_for_master_students/enya_and_johanna/genome_filepaths.tsv" 
file_paths_100k = [] # empty list
with open(filepaths, newline='', encoding='utf-8') as tsv_file:
    reader = csv.reader(tsv_file, delimiter='\t') # reads the file row by row
    for i, row in enumerate(reader): # we get i=line number, row=values in the row as a list
        if i == 100000:
            break
        file_paths_100k.append(row[0]) # takes column 0 from each row (we only have one column in the file)
'''    

# Bacterial genomes:
# get the first 100k lines in the file
file_paths_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_filepaths.tsv"
df = pd.read_csv(file_paths_file, sep="\t", header=None) 
file_paths = df[0].tolist()[1500_000:]

#print(file_paths_10)

# Open bacterial genome zip-files
# create one fasta file for all 100k genomes
output_fasta = "/storage/enyaa/REVISED/BLAST/FASTA/genomes_16.fasta"

with open(output_fasta, "w") as fasta_out:
    
    for path in file_paths: # goes through paths to the bacterial genome zip-files
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
blast_path = "/storage/enyaa/REVISED/BLAST/DB_16/blast_16" # name for BLAST database
subprocess.run(["makeblastdb", "-in", output_fasta, "-dbtype", "nucl", "-out", blast_path]) # runs this command in the terminal to create a database

end_time = time.time()

total_time = (end_time - start_time)/60

print(f"Database created with elapsed time: {total_time} minutes")


