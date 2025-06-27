import subprocess
import pandas as pd
import re
import tempfile
import gzip
import shutil
import concurrent.futures
import os
from tqdm import tqdm
'''
file = "/storage/enyaa/FINAL/filtered_bacteria.csv"
bacteria_df = pd.read_csv(file)
bacteria_id_set = set(bacteria_df['Bacteria_ID'])

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

        output_file = f"/storage/jolunds/FINAL/tRNA_SCAN/{bacteria_id}_trnascan.txt"

        # Run tRNAscan-SE
        command = ["tRNAscan-SE", "-B", "-o", output_file, tmp_fasta_path]
        subprocess.run(command, check=True)

    except subprocess.CalledProcessError as e:
        print(f"[!] {bacteria_id}: tRNAscan-SE failed: {e}")
    except Exception as e:
        print(f"[!] {bacteria_id}: Unexpected error: {e}")

num_parallel_jobs = 24
# Launch jobs in parallel
with concurrent.futures.ProcessPoolExecutor(max_workers=num_parallel_jobs) as executor:
    futures = [executor.submit(run_trnascan_job, bacteria_id, path) for bacteria_id, path in filtered_paths]
    
    for future in concurrent.futures.as_completed(futures):
        future.result()  # This will raise exceptions from the worker if any
'''


# Remove genomes with few tRNAs found
directory = "/storage/jolunds/FINAL/tRNA_SCAN/"
row_threshold = 35

# Loop through files, count rows and remove files with rows < 35
for filename in tqdm(os.listdir(directory), desc="Processing files"):
    filepath = os.path.join(directory, filename)
    if os.path.isfile(filepath):
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                for i, _ in enumerate(f, start=1):
                    if i >= row_threshold:
                        break
            if i < row_threshold:
                os.remove(filepath)
        except Exception as e:
            print(f"Error with {filename}: {e}")
                       