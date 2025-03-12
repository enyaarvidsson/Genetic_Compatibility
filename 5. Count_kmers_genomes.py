# Count k-mers for genomes
import subprocess
import os
import pandas as pd
import concurrent.futures
import time


# !!!!!!!!
# kör create dictionary innnan ny kmer counts
# ta bort genome_temp_remove och genome_kmer 
# innan jag startar ny
# !!!!!!!!


file_paths_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_filepaths.tsv"
df = pd.read_csv(file_paths_file, sep="\t", header=None) 
file_paths = df[0].tolist()[1000_000:1100_000] # extracts first column in df, converts to list, takes the first 10k filepaths

temp_dir = "/storage/enyaa/REVISED/KMER/genome_temp_remove" # Temporary directory
output_dir = "/storage/enyaa/REVISED/KMER/genome_kmer"

k = 5 # Decide k for k-mers
threads = 1 # Number for threads
num_parallel_jobs = 12 # Number of parallel jobs running, tror vi kan ha fler

start_time = time.time() # Starting time

running_jobs = [] 
with concurrent.futures.ProcessPoolExecutor(max_workers=num_parallel_jobs) as executor: # creates a pool of worker processes that can run in parallel
    futures = [] 
    
    for genome_fasta_file in file_paths:
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