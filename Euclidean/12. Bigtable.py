# Create a table with average compatiblity for each phylum
import pandas as pd
import time
import matplotlib.pyplot as plt
import os
import numpy as np
import seaborn as sns
from Bio import SeqIO

start_time = time.time()

# Load gene names
with open("/storage/jolunds/REVISED/gene_names.txt", "r") as f:
    all_genes = [line.strip() for line in f]

# Load full lineage 
path = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_lineage_df = pd.read_csv(path, sep="\t", header=None)
full_lineage_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
phylum_mapping = full_lineage_df[["Bacteria_ID", "Phylum", "Species"]] # only the bacteria_id and the respective phylum & species

# Loop through genes in all_genes 
big_table_list = []
table_genes_list = []   

for gene_name in all_genes: #sorted(all_genes): 
    
    if "/" in gene_name: # Just for look-up
        gene_name = gene_name.replace("/", "?")
    
    # Load euclidean df
    path = f"/storage/enyaa/REVISED/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl"
    gene_euclidean_df = pd.read_pickle(path).T.reset_index() # Switch to long format 
    gene_euclidean_df.columns = ['Bacteria_ID', 'Euclidean_distance']
 
    # Add phylum column
    phylum_euclidean_df = gene_euclidean_df.merge(phylum_mapping, on=['Bacteria_ID'], how='inner')
    #print(phylum_euclidean_df.head())

    # Find 6 top phylum & filter for top phyla
    top_phyla = phylum_euclidean_df["Phylum"].value_counts().head(6)
    phylum_euclidean_df = phylum_euclidean_df[phylum_euclidean_df['Phylum'].isin(top_phyla.index)]
    
    # Calculate mean and standard deviation for each phylum and gene
    phylum_stats = phylum_euclidean_df.groupby(["Phylum"])["Euclidean_distance"].agg(['mean', 'std']).fillna(0)
    phylum_stats_reset = phylum_stats.reset_index()
        # has the genes as rows, and shows mean and std for each phylum

    # Load taxonomy results and count matches
    taxonomy_file = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/taxonomy_results_{gene_name}.csv"
    if os.path.exists(taxonomy_file):
        taxonomy_df = pd.read_csv(taxonomy_file)
        taxonomy_df.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        taxonomy_df = taxonomy_df.drop_duplicates() # Remove duplicates (when we have multiple matches in one genome)
        
        matching_df = taxonomy_df[['Bacteria_ID']]
        phylum_euclidean_df = phylum_euclidean_df.merge(matching_df.assign(Match_status="Match"), on="Bacteria_ID", how="left").fillna("No_match") # here a new column is added to euclidean_gene_df called "Match_status" and it says if there are Match
        phylum_counts = taxonomy_df.groupby('Phylum').size().reset_index(name="Num_matches") # Count num matches in each phylum
    else: # If taxonomy file do not exist
        phylum_counts = pd.DataFrame({"Phylum": top_phyla.index.tolist(), "Num_matches": [0] * len(top_phyla)})
        phylum_euclidean_df["Match_status"] = "No_match"
    
    # Merge stats and num matches
    gene_table_df = pd.merge(phylum_stats_reset, phylum_counts, on=["Phylum"], how="left").fillna(0) # Add num matches 
    
    # Total mean, minimum value, minimum value for match + species, maximum value for match + species
    total_mean = phylum_euclidean_df['Euclidean_distance'].mean()
    minimum_eu = phylum_euclidean_df['Euclidean_distance'].min()
    maximum_eu = phylum_euclidean_df['Euclidean_distance'].max()
    
    match_df = phylum_euclidean_df[phylum_euclidean_df['Match_status'] == 'Match']
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
        'Total_eu_mean': [total_mean],
        'Min_eu': [minimum_eu],
        'Max_eu': [maximum_eu],
        'Min_eu_match': [min_match],
        'Min_Species': [min_species],
        'Max_eu_match': [max_match],
        'Max_Species': [max_species],
        'Mean_eu_matches': [mean_match],
        'Num_phyla_match': [unique_phyla_count]
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
    gene_table_df['Phylum'] = pd.Categorical(gene_table_df['Phylum'], categories=top_phyla.index, ordered=True) 
    gene_table_df = gene_table_df.sort_values('Phylum') # Sort phylum from largest to smallest
    
    # Append results
    big_table_list.append(gene_table_df)
    table_genes_list.append(table_gene_df)

# Concat
big_table_df = pd.concat(big_table_list).reset_index(drop=True)
print(big_table_df)
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
print(table_genes_df)


# Save
path_bigtable = "/storage/jolunds/REVISED/KMER/big_table.csv"
big_table_df.to_csv(path_bigtable, index=False)
path_table_genes = "/storage/jolunds/REVISED/KMER/table_genes.csv"
table_genes_df.to_csv(path_table_genes, index=False)

end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Done creating big table with elapsed time: {total_time} minutes")

'''
# Create scatterplot
start_time = time.time()

big_table_df = pd.read_csv("/storage/jolunds/REVISED/KMER/big_table.csv")

big_table_df["Num_matches"] = big_table_df["Num_matches"].clip(upper=5000)

# Använda facet grid för att göra en figur med en plot för varje phylum
g = sns.FacetGrid(big_table_df, col="Phylum", col_wrap = 3)
g.map_dataframe(sns.scatterplot, x="Num_matches", y="Mean")
g.set_axis_labels("Number of matches", "Mean euclidean distance",)

# Save plot
plt.savefig("/home/jolunds/newtest/scatterplot_big_table_5.png")
plt.close(g.figure)
top2_per_phylum = big_table_df.groupby("Phylum").apply(lambda x: x.nlargest(2, "Mean")).reset_index(drop=True)

# Display the results
print(top2_per_phylum[["Phylum", "Gene_name", "Mean"]])
end_time = time.time()
total_time = (end_time - start_time)/60


print(f"Done creating big table with elapsed time: {total_time} minutes")
'''


