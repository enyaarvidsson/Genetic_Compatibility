import pandas as pd
import time
from tqdm import tqdm

start_time = time.time()
# Get the gene IDs and bacteria IDs from the BLAST results
with open(f"/storage/jolunds/FINAL/BLAST/filtered_blast_results.tsv") as f: # opens the file for reading, stores in f
    lines = f.readlines()[1:] # reads the lines in the file and stores them in a list "lines"

gene_header = [line.split("\t")[0] for line in lines]

gene_names = [line.split("|")[5] for line in gene_header] # Takes gene names

bacteria_id_seq = [line.split("\t")[1] for line in lines] # Takes the ID for the sequences in the bacterial genomes
bacteria_id = [line.split("|")[0] for line in bacteria_id_seq] # Removes sequence ID, so we just have bacteria IDs

blast_results_df = pd.DataFrame({
        "Gene_name": gene_names,
        "Bacteria_ID": bacteria_id
    })

# Load gene names
with open("/storage/enyaa/FINAL/gene_names.txt", "r") as f:
    all_genes = [line.strip() for line in f]

# Read taxonomy file
taxonomy_file = "/storage/enyaa/FINAL/filtered_bacteria.csv"
taxonomy_df = pd.read_csv(taxonomy_file) 

gene_count = 0
for gene_name in tqdm(all_genes, desc="Processing genes"):
    gene_count += 1

    if gene_name in blast_results_df["Gene_name"].values:
        matched_bacteria = blast_results_df[blast_results_df["Gene_name"] == gene_name]["Bacteria_ID"]
    else:
        matched_bacteria = []
    
    gene_df = taxonomy_df.copy()
    gene_df["Match_status"] = gene_df["Bacteria_ID"].apply(
        lambda x: "Match" if x in set(matched_bacteria) else "Non-match"
    )

    if "/" in gene_name: # Just for look-up
        gene_name = gene_name.replace("/", "?")
    gene_df.to_csv(f"/storage/jolunds/FINAL/MATCH_INFO/{gene_name}_match_info.tsv", sep="\t", index=False)

end_time = time.time()
total_time = (end_time -start_time)/60

print(f"Done creating match info for {gene_count} genes, with elapsed time {total_time} minutes")
