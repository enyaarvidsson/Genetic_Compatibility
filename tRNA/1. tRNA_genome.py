import subprocess
import pandas as pd
import re
import tempfile
import gzip
import shutil
import os
import time
import concurrent.futures

start_time = time.time()

file = "/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_tet(Q).pkl"
downsampled_df = pd.read_pickle(file)

bacteria_id_set = set(downsampled_df['Bacteria_ID'])

pattern = re.compile(r'(GCA_\d+\.\d+)')

'''filtered_paths = []
filepaths = "/storage/shared/data_for_master_students/enya_and_johanna/genome_filepaths.tsv"
with open(filepaths) as file:
    
    for line in file:
        path = line.strip()
        match = pattern.search(path)
        if match:
            bacteria_id = match.group(1)
            if bacteria_id in bacteria_id_set:
                filtered_paths.append((bacteria_id, path))
'''
def run_trnascan_job(bacteria_id, path):
    tmp_fasta_path = None
    try:
        # Create a temporary unzipped FASTA file
        with tempfile.NamedTemporaryFile(suffix=".fna", delete=False, dir='/storage/jolunds/temp_dir') as tmp_fasta:
            with gzip.open(path, 'rb') as f_in:
                shutil.copyfileobj(f_in, tmp_fasta)
            tmp_fasta_path = tmp_fasta.name

        output_file = f"/storage/jolunds/REVISED/tRNA/{bacteria_id}_trnascan.txt"

        # Run tRNAscan-SE
        command = ["tRNAscan-SE", "-B", "-o", output_file, tmp_fasta_path]
        subprocess.run(command, check=True)

    except subprocess.CalledProcessError as e:
        print(f"[!] {bacteria_id}: tRNAscan-SE failed: {e}")
    except Exception as e:
        print(f"[!] {bacteria_id}: Unexpected error: {e}")


'''num_parallel_jobs = 24
# Launch jobs in parallel
with concurrent.futures.ProcessPoolExecutor(max_workers=num_parallel_jobs) as executor:
    futures = [executor.submit(run_trnascan_job, bacteria_id, path) for bacteria_id, path in filtered_paths]
    
    for future in concurrent.futures.as_completed(futures):
        future.result()  # This will raise exceptions from the worker if any
'''

# Koden gav error-meddelanden så många filer är tomma
# Hitta bacteria ids som ej har gjorts trna scan på

# Get all filenames in the directory
files_in_dir = set(os.listdir("/storage/jolunds/REVISED/tRNA/"))

# Extract the GCA IDs from the filenames (assuming they start like GCA_12345.1_...)
completed_bacteria_ids = set()
for filename in files_in_dir:
    if filename.startswith("GCA_") and filename.endswith(".txt"):
        parts = filename.split("_")
        if len(parts) >= 2:
            bacteria_id = f"{parts[0]}_{parts[1]}"
            completed_bacteria_ids.add(bacteria_id)

# Find bacteria IDs that are NOT present in the directory
error_bacteria_ids = bacteria_id_set - completed_bacteria_ids

pattern = re.compile(r'(GCA_\d+\.\d+)')
error_file_paths = []
filepaths = "/storage/shared/data_for_master_students/enya_and_johanna/genome_filepaths.tsv"
with open(filepaths) as file:
    
    for line in file:
        path = line.strip()
        match = pattern.search(path)
        if match:
            bacteria_id = match.group(1)
            if bacteria_id in error_bacteria_ids:
                error_file_paths.append((bacteria_id, path))

num_parallel_jobs = 2
# Launch jobs in parallel
with concurrent.futures.ProcessPoolExecutor(max_workers=num_parallel_jobs) as executor:
    futures = [executor.submit(run_trnascan_job, bacteria_id, path) for bacteria_id, path in error_file_paths]
    
    for future in concurrent.futures.as_completed(futures):
        future.result()  # This will raise exceptions from the worker if any

end_time = time.time()
total_time = (end_time-start_time)/60
print(f"Done runnning tRNAscan, total time {total_time} minutes")