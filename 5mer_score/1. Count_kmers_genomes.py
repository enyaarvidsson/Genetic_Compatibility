# Count k-mers (5-mers) for the bacterial genomes

import subprocess
import os
import pandas as pd
import concurrent.futures
import time
from tqdm import tqdm


# !!!!!!!!
# kör create dictionary innnan ny kmer counts
# ta bort genome_temp_remove och genome_kmer 
# innan jag startar ny
# !!!!!!!!


# Filtered bacteria -----------------
filtered_bacteria_df = pd.read_csv("/storage/enyaa/FINAL/filtered_bacteria.csv")
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



temp_dir = "/storage/enyaa/FINAL/KMER/genome_temp_remove" # Temporary directory
output_dir = "/storage/enyaa/FINAL/KMER/genome_kmer"

k = 5 # Decide k for k-mers
threads = 1 # Number of threads
num_parallel_jobs = 12 # Number of parallel jobs running

start_time = time.time() 

running_jobs = [] 
with concurrent.futures.ProcessPoolExecutor(max_workers=num_parallel_jobs) as executor: # creates a pool of worker processes that can run in parallel
    futures = [] 
    
    for genome_fasta_file in tqdm(file_paths, desc="Processing genomes"):
        filename = os.path.basename(genome_fasta_file)
        bacteria_id = "_".join(filename.split("_")[:2])
    
        output_prefix = os.path.join(temp_dir, f"{bacteria_id}_kmc")
        output_txt = os.path.join(output_dir, f"{bacteria_id}_kmc_counts.txt")
        running_jobs.append((output_prefix, output_txt)) # saves output_prefix and output_txt in "running_jobs"
        
        futures.append(executor.submit(
            subprocess.run, 
            ["kmc", "-fm", f"-k{k}", f"-t{threads}", "-ci1", "-cs500000",
             genome_fasta_file, str(output_prefix), temp_dir], 
            check=True
        ))
    # Wait for all jobs to complete before dumping k-mer counts
    concurrent.futures.wait(futures)

    # Now, run kmc_dump in parallel
    dump_futures = [
        executor.submit(subprocess.run, ["kmc_dump", str(output_prefix), str(output_txt)], check=True)
        for output_prefix, output_txt in running_jobs
    ]

    concurrent.futures.wait(dump_futures)

end_time = time.time()

total_time = (end_time - start_time)/60

print(f"Done processing KMC for all genomes with elapsed time: {total_time} minutes")
