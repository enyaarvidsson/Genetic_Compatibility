import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

gene_name = "AAC6_30_AAC6_Ib"
tRNA_score = "tRNA_score_one_sided"

if "/" in gene_name:
    gene_name = gene_name.replace("/", "?")
    
file_path = f"/storage/jolunds/REVISED/tRNA/tRNA_score/tRNA_score_{gene_name}.csv"
tRNA_score_df = pd.read_csv(file_path)

# Load taxonomy results and count matches
taxonomy_file = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/taxonomy_results_{gene_name}.csv"
if os.path.exists(taxonomy_file):
    taxonomy_df = pd.read_csv(taxonomy_file)
    taxonomy_df.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    taxonomy_df = taxonomy_df.drop_duplicates() # Remove duplicates (when we have multiple matches in one genome)
    matching_df = taxonomy_df[['Bacteria_ID']]
    tRNA_score_df = tRNA_score_df.merge(matching_df.assign(Match_status="Match"), on="Bacteria_ID", how="left").fillna("No_match") # here a new column is added to euclidean_gene_df called "Match_status" and it says if there are Match
    matches = 1
else: # If taxonomy file do not exist
    tRNA_score_df["Match_status"] = "No_match"
    matches = 0
    
no_match_count = tRNA_score_df["Match_status"].value_counts().get("Match", 0)    
if no_match_count == 0:
        print("No matches for gene:", gene_name)
        matches = 0

# Filter for top phyla
top_phyla = tRNA_score_df["Phylum"].value_counts().head(6)
tRNA_score_df = tRNA_score_df[tRNA_score_df['Phylum'].isin(top_phyla.index)]

nr_bins = 30
min_value = tRNA_score_df[tRNA_score].min()
max_value = tRNA_score_df[tRNA_score].max() + 0.001 # so all values fall inside the max_value
bin_edges = np.linspace(min_value, max_value, nr_bins + 1)

# DOWNSAMPLE NO_MATCH
downsampled_no_matches = [] # will become a list of dataframes
tRNA_score_df = tRNA_score_df.copy()
tRNA_score_df["Phylum"] = tRNA_score_df["Phylum"].astype(str)

for phylum, phylum_df in tRNA_score_df.groupby("Phylum"): # phylum - name of phylum, phylum_df - df that has only rows from that phylum
    match_count = (phylum_df["Match_status"] == "Match").sum()
    # create a df for match and a df for no_match
    matches_phylum_df = phylum_df[phylum_df["Match_status"] == "Match"]
    no_matches_phylum_df = phylum_df[phylum_df["Match_status"] == "No_match"]  

    # How many no_matches to keep
    if matches == 1: # if matches exists
        if match_count == 0:
            keep_size = 10000
        elif match_count < 100:
            keep_size = 2000
        elif match_count < 3000:
            keep_size = 10000
        else:
            keep_size = match_count * 3
    else:
        keep_size = len(phylum_df) 

    # Downsample no_matches
    if keep_size > len(no_matches_phylum_df): # make sure np.random.choice works, so the keep_size is not bigger than the population size
        keep_size = len(no_matches_phylum_df)
    downsampled_no_matches_phylum_df = no_matches_phylum_df.iloc[np.random.choice(len(no_matches_phylum_df), keep_size, replace=False)]

    # Append both "Match" and downsampled "No_match" bacteria
    downsampled_no_matches.append(pd.concat([matches_phylum_df, downsampled_no_matches_phylum_df]))

# Combine all phyla into the final df
tRNA_downsampled_df = pd.concat(downsampled_no_matches, ignore_index=True)

phylum_counts = tRNA_downsampled_df['Phylum'].value_counts()  
phylum_counts = phylum_counts.reindex(top_phyla.index) # In same order as top_phyla

# matches in each phylum
matches_phylum_counts = tRNA_downsampled_df[tRNA_downsampled_df['Match_status'] == 'Match'].groupby('Phylum').size()
matches_phylum_counts = matches_phylum_counts.reindex(top_phyla.index).fillna(0).astype(int)


# HISTOGRAM
g = sns.FacetGrid(tRNA_downsampled_df, col="Phylum", col_order=top_phyla.index, sharey=False, col_wrap=3, height=4, aspect=1.2)
g.map_dataframe(sns.histplot, x=tRNA_score, hue = "Match_status", hue_order=["No_match", "Match"], multiple="stack", bins=bin_edges)
g.set_axis_labels("tRNA score", "Number of Bacteria")

for ax, phylum in zip(g.axes.flat, phylum_counts.index):
    ax.set_title(f"{phylum} (n={phylum_counts[phylum]}, m={matches_phylum_counts[phylum]})")
        
g.set(xlim=(min_value - 0.01, max_value + 0.01))
    
plt.subplots_adjust(top=0.85)

 # Add title
if "?" in gene_name:
    gene_name = gene_name.replace("?", "/")
    
if matches == 1: # if matches exists
    plt.figtext(0.5, 0.95, f"Gene name: {gene_name}", ha="center", fontsize=14)
else:
    plt.figtext(0.5, 0.95, f"Gene name: {gene_name} - NO MATCHES", ha="center", fontsize=14)     

plt.savefig(f'/home/enyaa/gene_genome/histogram_{gene_name}.png')     
plt.close(g.figure)

print(f"Histogram created for {gene_name}!")