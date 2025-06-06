# Run BLAST between 6k genes and 70k bacteria - to find matches

import subprocess


blast_db = f"/storage/enyaa/FINAL/BLAST/DB/blast"
query_file = "/storage/enyaa/nucleotide_fasta_protein_homolog_model.fasta" # all 6k genes
store_results = f"/storage/enyaa/FINAL/BLAST/blast_results.txt"
subprocess.run(["blastn", "-query", query_file, "-db", blast_db, "-out", store_results, "-perc_identity", "95", 
                "-max_target_seqs", "100000", "-num_threads", "8", "-evalue", "1e-5", "-qcov_hsp_perc", "90", "-best_hit_score_edge", "0.1", "-outfmt", "6"],
                check=True) # "-outfmt", "6" - means that output is in a tabular format (tab-separated file)    
print(f"Blast_results created!")

