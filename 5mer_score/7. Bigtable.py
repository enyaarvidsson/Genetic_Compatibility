# Create a table with average compatiblity for each phylum
import pandas as pd
import time
import os
from Bio import SeqIO
from tqdm import tqdm

start_time = time.time()

# Load gene names
with open("/storage/enyaa/FINAL/gene_names.txt", "r") as f:
    all_genes = [line.strip() for line in f]

# Loop through genes in all_genes 
big_table_list = []
table_genes_list = []   

for gene_name in tqdm(all_genes, desc="Processing genes"): #sorted(all_genes): 
    
    if "/" in gene_name: # Just for look-up
        gene_name = gene_name.replace("/", "?")
    
    # Load 5mer score
    path = f"/storage/enyaa/FINAL/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl" #original
    
    if not os.path.exists(path):
        continue 
    
    gene_euclidean_df = pd.read_pickle(path)
    
    # Calculate mean and standard deviation for each phylum and gene
    phylum_stats = gene_euclidean_df.groupby(["Phylum"])["Euclidean_distance"].agg(['mean', 'std']).fillna(0)
    phylum_stats_reset = phylum_stats.reset_index()

    # Calculate number of matches
    phylum_counts = gene_euclidean_df[gene_euclidean_df["Match_status"] == "Match"].groupby("Phylum").size().reset_index(name="Num_matches")
    
    # Merge stats and num matches
    gene_table_df = pd.merge(phylum_stats_reset, phylum_counts, on=["Phylum"], how="left").fillna(0) # Add num matches 
    
    # Total mean, minimum value, minimum value for match + species, maximum value for match + species
    total_mean = gene_euclidean_df['Euclidean_distance'].mean()
    minimum_eu = gene_euclidean_df['Euclidean_distance'].min()
    maximum_eu = gene_euclidean_df['Euclidean_distance'].max()
    
    match_df = gene_euclidean_df[gene_euclidean_df['Match_status'] == 'Match']
    
    if match_df.empty:
        # Per gene
        min_match = "No matches"
        min_species = "No matches"
        max_match = "No matches"
        max_species = "No matches"
        mean_match = "No matches"
        unique_phyla_count = 0
        
        # Per phylum
        phylum_min_match = pd.Series(dtype=float)
        phylum_min_species = pd.Series(dtype=str)
        phylum_max_match = pd.Series(dtype=float)
        phylum_max_species = pd.Series(dtype=str)
        phylum_mean_match = pd.Series(dtype=float)
        
    else:   
        # Per gene
        min_row = match_df.loc[match_df['Euclidean_distance'].idxmin(), ['Euclidean_distance', 'Species']]
        min_match = min_row['Euclidean_distance']
        min_species = min_row['Species']
        
        max_row = match_df.loc[match_df['Euclidean_distance'].idxmax(), ['Euclidean_distance', 'Species']]
        max_match = max_row['Euclidean_distance']
        max_species = max_row['Species']
        
        mean_match = match_df['Euclidean_distance'].mean()
        unique_phyla_count = match_df['Phylum'].nunique()
        
        # Per phylum
        phylum_min_match = match_df.groupby('Phylum')['Euclidean_distance'].min()
        phylum_max_match = match_df.groupby('Phylum')['Euclidean_distance'].max()
        phylum_mean_match = match_df.groupby('Phylum')['Euclidean_distance'].mean()
    
        
        phylum_min_species = match_df.loc[match_df.groupby('Phylum')['Euclidean_distance'].idxmin(), ['Phylum', 'Species']].set_index('Phylum')['Species']
        phylum_max_species = match_df.loc[match_df.groupby('Phylum')['Euclidean_distance'].idxmax(), ['Phylum', 'Species']].set_index('Phylum')['Species']
    
    # Create tables
    if "?" in gene_name:
        gene_name = gene_name.replace("?", "/")
    table_gene_df = pd.DataFrame({
        'Gene_name': [gene_name],
        'Mean': [total_mean],
        'Min': [minimum_eu],
        'Max': [maximum_eu],
        'Min_match': [min_match],
        'Min_species': [min_species],
        'Max_eu_match': [max_match],
        'Max_species': [max_species],
        'Mean_matches': [mean_match],
        'Num_phyla': [unique_phyla_count]
    })
    
    gene_table_df['Min_match'] = gene_table_df['Phylum'].map(phylum_min_match).fillna("No matches")
    gene_table_df['Min_species'] = gene_table_df['Phylum'].map(phylum_min_species).fillna("No matches")
    gene_table_df['Max_match'] = gene_table_df['Phylum'].map(phylum_max_match).fillna("No matches")
    gene_table_df['Max_species'] = gene_table_df['Phylum'].map(phylum_max_species).fillna("No matches")
    gene_table_df['Mean_match'] = gene_table_df['Phylum'].map(phylum_mean_match).fillna("No matches")
    
    # Changes for the eye    
    gene_table_df.insert(0, 'Gene_name', gene_name) # Add gene name column
    gene_table_df.rename(columns={'mean': 'Mean', 'std': 'Std'}, inplace=True) # Capitalise mean and std
    gene_table_df['Num_matches'] = gene_table_df['Num_matches'].astype(int) # Make num matches integers
    top_phyla = ["Pseudomonadota", "Bacillota", "Bacteroidota", "Cyanobacteriota", "Actinomycetota", "Campylobacterota"]
    gene_table_df['Phylum'] = pd.Categorical(gene_table_df['Phylum'], categories=top_phyla, ordered=True) 
    gene_table_df = gene_table_df.sort_values('Phylum') # Sort phylum from largest to smallest
    
    # Append results
    big_table_list.append(gene_table_df)
    table_genes_list.append(table_gene_df)

# Concat
big_table_df = pd.concat(big_table_list).reset_index(drop=True)
table_genes_df = pd.concat(table_genes_list).reset_index(drop=True)

# Gene length
filepath = "/storage/enyaa/nucleotide_fasta_protein_homolog_model.fasta"

gene_lengths = {}
for record in SeqIO.parse(filepath, "fasta"):
    header = record.description
    
    positions = header.split("|")[3] 
    start = int(positions.split("-")[0])
    end = int(positions.split("-")[1])
    gene_length = abs(end-start) + 1
    gene_name = header.split("|")[-1].split(" ")[0]
    
    gene_lengths[gene_name] = gene_length
    
gene_lengths_df = pd.DataFrame(list(gene_lengths.items()), columns=['Gene_name', 'Gene_length'])

# Add gene length
table_genes_df = table_genes_df.merge(gene_lengths_df, on='Gene_name', how='left')

# Save
path_bigtable = "/storage/jolunds/FINAL/5mer_score_bigtable.csv"
big_table_df.to_csv(path_bigtable, index=False)
path_table_genes = "/storage/jolunds/FINAL/5mer_score_table_genes.csv"
table_genes_df.to_csv(path_table_genes, index=False)

end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Done creating big tables with elapsed time: {total_time} minutes")
