# Run BLAST

#from Bio import SeqIO # read and write sequence files
import subprocess
import time

'''
# For taking out 10 genes: 

# Filepath for FASTA file with the genes
file_path = "/storage/enyaa/nucleotide_fasta_protein_homolog_model.fasta"

# take out 10 genes
records = SeqIO.parse(file_path, "fasta")
genes10 = [next(records) for _ in range(10)]

query_file = "/home/enyaa/gene_genome/genes10.fasta"
with open(query_file, "w") as out_file:
    SeqIO.write(genes10, out_file, "fasta")
'''

start_time = time.time()

# Run BLASTn between 6000 genes and 100k genomes
blast_db = "/storage/enyaa/REVISED/BLAST/DB_15/blast_15"
query_file = "/storage/enyaa/nucleotide_fasta_protein_homolog_model.fasta" # all 6k genes
store_results = "/storage/enyaa/REVISED/BLAST/BLAST_RESULTS/blast_results_15.txt"
subprocess.run(["blastn", "-query", query_file, "-db", blast_db, "-out", store_results, "-perc_identity", "95", 
                "-max_target_seqs", "100000", "-num_threads", "8", "-evalue", "1e-5", "-outfmt", "6"],
                check=True) # "-outfmt", "6" - means that output is in a tabular format (tab-separated file)    

end_time = time.time()
total_time = (end_time -start_time)/60

print(f"Done running BLASTn with elapsed time {total_time} minutes")

'''
# To check if we got any matches:

with open("blast_results_10k.txt", "r") as f: # opens the BLAST search file for reading, with means it closes automatically
    results = f.read() # stores the reading in results

# Print results if there exists any matches
if not results:
    print("No significant matches found.")
else:
    print("BLAST search completed. Results saved in blast_results_10k.txt.")
'''