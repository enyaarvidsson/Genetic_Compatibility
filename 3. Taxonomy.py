# 3. Taxonomy

import pandas as pd

# Get the gene IDs and bacteria IDs from the BLAST results
with open("/storage/enyaa/REVISED/BLAST/BLAST_RESULTS/blast_results_9.txt") as f: # opens the file for reading, stores in f
    lines = f.readlines() # reads the lines in the file and stores them in a list "lines"

gene_header = [line.split("\t")[0] for line in lines]
gene_name = [line.split("|")[5] for line in gene_header] # Takes gene names
#print(gene_name[:10])

bacteria_id_seq = [line.split("\t")[1] for line in lines] # Takes the ID for the sequences in the bacterial genomes
bacteria_id = [line.split("|")[0] for line in bacteria_id_seq] # Removes sequence ID, so we just have bacteria IDs

# Read taxonomy file
taxonomy_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
taxonomy_df = pd.read_csv(taxonomy_file, sep="\t", header=None) 

# Find matching genomes by bacteria ID
matching_rows = taxonomy_df[taxonomy_df.iloc[:, 0].isin(bacteria_id)] # Matches from taxonomy and results, gets the whole row in taxonomy
    # taxonomy_df.iloc[:, 0] extracts the first column as a Pandas Series.
    # .isin(bacteria_id) checks for each value in the first column whether it exists in bacteria_id, 
        # returning a Boolean Series (True or False)
    # taxonomy_df[...] selects only the rows where the condition inside is True.

matching_rows.insert(0, 'Gene_ID', gene_name[:len(matching_rows)]) # Inserts a column, in position 0, with the gene ID in matching rows
    # gene_id[:len(matching_rows)] slices the first len(matching_rows) elements from gene_id, 
    # ensuring that the number of inserted values matches the number of rows in matching_rows
    
# Store results
store_results = "/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_9.csv"
matching_rows.to_csv(store_results, index=False, header=False) # Creates file with the results



