# Creates a histogram for a specific gene, for the tRNA score

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np


gene_name = "tet(Q)" 


np.random.seed(42)

if "/" in gene_name:
    gene_name = gene_name.replace("/", "?")

file_path = f"./FINAL/tRNA_SCORE/{gene_name}_tRNA_score.csv"
tRNA_score_df = pd.read_csv(file_path)
    # 67098 rows (bacteria), with columns: Bacteria_ID, tRNA_score, Worst_case_tRNA, Worst_codon, taxonomic info (multiple columns), Match_status

matches = 1    
no_match_count = tRNA_score_df["Match_status"].value_counts().get("Match", 0)    
if no_match_count == 0:
    print("No matches for gene:", gene_name)
    matches = 0

# Filter for top phyla
top_phyla = tRNA_score_df["Phylum"].value_counts().head(6)

# Same order as 5mer score
phyla_order = ["Pseudomonadota", "Bacillota", "Actinomycetota", "Bacteroidota", "Cyanobacteriota", "Campylobacterota"]
top_phyla = top_phyla.loc[phyla_order]

tRNA_score_df = tRNA_score_df[tRNA_score_df['Phylum'].isin(top_phyla.index)]


# DOWNSAMPLE NON-MATCH
downsampled_no_matches = [] # will become a list of dataframes
tRNA_score_df = tRNA_score_df.copy()
tRNA_score_df["Phylum"] = tRNA_score_df["Phylum"].astype(str)

for phylum, phylum_df in tRNA_score_df.groupby("Phylum"): # phylum - name of phylum, phylum_df - df that has only rows from that phylum
    match_count = (phylum_df["Match_status"] == "Match").sum()
    # create a df for match and a df for no_match
    matches_phylum_df = phylum_df[phylum_df["Match_status"] == "Match"]
    no_matches_phylum_df = phylum_df[phylum_df["Match_status"] == "Non-match"]  

    # How many no_matches to keep
    if matches == 1: # if matches exists
        if match_count == 0:
            keep_size = 5000
        elif match_count < 100:
            keep_size = 1000
        elif match_count < 3000:
            keep_size = 5000
        else:
            keep_size = match_count * 3
    else:
        keep_size = len(phylum_df) 

    # Downsample non-matches
    if keep_size > len(no_matches_phylum_df): # make sure np.random.choice works, so the keep_size is not bigger than the population size
        keep_size = len(no_matches_phylum_df)
    downsampled_no_matches_phylum_df = no_matches_phylum_df.iloc[np.random.choice(len(no_matches_phylum_df), keep_size, replace=False)]

    # Append both "Match" and downsampled "Non-match" bacteria
    downsampled_no_matches.append(pd.concat([matches_phylum_df, downsampled_no_matches_phylum_df]))

# Combine all phyla into the final df
tRNA_downsampled_df = pd.concat(downsampled_no_matches, ignore_index=True)

phylum_counts = tRNA_downsampled_df['Phylum'].value_counts()  
phylum_counts = phylum_counts.reindex(top_phyla.index) # In same order as top_phyla

# matches in each phylum
matches_phylum_counts = tRNA_downsampled_df[tRNA_downsampled_df['Match_status'] == 'Match'].groupby('Phylum').size()
matches_phylum_counts = matches_phylum_counts.reindex(top_phyla.index).fillna(0).astype(int)

nr_bins = 30
min_value = tRNA_downsampled_df["tRNA_score"].min()
max_value = tRNA_downsampled_df["tRNA_score"].max() 
bin_edges = np.linspace(min_value, max_value, nr_bins + 1)


# HISTOGRAM
g = sns.FacetGrid(tRNA_downsampled_df, col="Phylum", col_order=top_phyla.index, sharey=False, col_wrap=3, height=4, aspect=1.2)
g.map_dataframe(sns.histplot, x="tRNA_score", hue = "Match_status", hue_order=["Non-match", "Match"], multiple="stack", bins=bin_edges)
g.set_axis_labels("tRNA score", "") # Number of bacteria

for ax, phylum in zip(g.axes.flat, phylum_counts.index):
    ax.set_title(f"{phylum} (n={phylum_counts[phylum]}, m={matches_phylum_counts[phylum]})", fontsize=15)
    #ax.set_title(f"{phylum}", fontsize=15)
    ax.set_xlabel("tRNA score", fontsize=15)
    ax.tick_params(axis='both', labelsize=13)    

g.set(xlim=(min_value-0.003, max_value))
    
plt.subplots_adjust(top=0.85)

# Add title
if "?" in gene_name:
    gene_name = gene_name.replace("?", "/")

if matches == 1: # if matches exists
    plt.figtext(0.5, 0.95, f"Gene name: {gene_name}", ha="center", fontsize=14)
else:
    plt.figtext(0.5, 0.95, f"Gene name: {gene_name} - NO MATCHES", ha="center", fontsize=14)     

#plt.tight_layout() # in report
plt.savefig(f'./FINAL/histogram_tRNA_score_{gene_name}.png')     
plt.close(g.figure)

print(f"tRNA score histogram created for {gene_name}")
