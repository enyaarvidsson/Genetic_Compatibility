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

filtered_paths = []
filepaths = "/storage/shared/data_for_master_students/enya_and_johanna/genome_filepaths.tsv"
with open(filepaths) as file:
    
    for line in file:
        path = line.strip()
        match = pattern.search(path)
        if match:
            bacteria_id = match.group(1)
            if bacteria_id in bacteria_id_set:
                filtered_paths.append((bacteria_id, path))

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
    finally:
        # Clean up the temporary file
        if tmp_fasta_path and os.path.exists(tmp_fasta_path):
            os.remove(tmp_fasta_path)

num_parallel_jobs = 24
# Launch jobs in parallel
with concurrent.futures.ProcessPoolExecutor(max_workers=num_parallel_jobs) as executor:
    futures = [executor.submit(run_trnascan_job, bacteria_id, path) for bacteria_id, path in filtered_paths]
    
    for future in concurrent.futures.as_completed(futures):
        future.result()  # This will raise exceptions from the worker if any

end_time = time.time()
total_time = (end_time-start_time)/60
print(f"Done runnning tRNAscan, total time {total_time} minutes")

