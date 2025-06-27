
import matplotlib.pyplot as plt
import pandas as pd
import os
import seaborn as sns
from scipy.stats import spearmanr
import numpy as np


gene_name = "NDM-1" 
np.random.seed(42)

if "/" in gene_name:
    gene_name = gene_name.replace("/", "?")


# tRNA SCORE ------------------
filepath_tRNA = f"/storage/jolunds/FINAL/tRNA_SCORE/{gene_name}_tRNA_score.csv"
tRNA_score_df = pd.read_csv(filepath_tRNA)
    
no_match_count = tRNA_score_df["Match_status"].value_counts().get("Match", 0)    
if no_match_count == 0:
    print("No matches for gene:", gene_name)
    matches = 0

# ONLY MATHCES ----------------
#tRNA_score_df = tRNA_score_df[tRNA_score_df["Match_status"] == 'Match']


# DOWNSAMPLE NO-MATCHES -------
tRNA_score_downsampled = [] # will become a list of dataframes

match_count = (tRNA_score_df["Match_status"] == "Match").sum()

# create a df for match and a df for no_match
matches_df = tRNA_score_df[tRNA_score_df["Match_status"] == "Match"]
no_matches_df = tRNA_score_df[tRNA_score_df["Match_status"] == "Non_match"]  

# How many no_matches to keep
if matches == 1: # if matches exists
    if match_count < 100:
        keep_size = 100
    else:
        keep_size = match_count
else:
    print("No matches") 

# Downsample no_matches
if keep_size > len(no_matches_df): 
    keep_size = len(no_matches_df)
downsampled_no_matches_df = no_matches_df.iloc[np.random.choice(len(no_matches_df), keep_size, replace=False)]

# Append both "Match" and downsampled "No_match" bacteria
tRNA_score_downsampled.append(pd.concat([matches_df, downsampled_no_matches_df]))

tRNA_score_downsampled_df = pd.concat(tRNA_score_downsampled, ignore_index=True)


# EUCLIDEAN DISTANCE ----------
filepath_eu = f"/storage/enyaa/FINAL/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl"
euclidean_df = pd.read_pickle(filepath_eu)

# MERGE THEM TOGETHER ----------
tRNA_and_euclidean_df = pd.merge(tRNA_score_downsampled_df, euclidean_df, on='Bacteria_ID', how='inner')


# SPEARMAN CORRELATION --------
correlation, p_value = spearmanr(tRNA_and_euclidean_df['Euclidean_distance'], tRNA_and_euclidean_df['tRNA_score'])

print(f"Spearman correlation: {correlation}")
print(f"P-value: {p_value}")


# SCATTERPLOT -----------
sns.scatterplot(
    data=tRNA_and_euclidean_df,
    x='Euclidean_distance',
    y='tRNA_score',   
    hue='Match_status',
    hue_order=["Non_match", "Match"],
    alpha=1,
    s=10
    #color='darkorange'
)

plt.xlabel('Euclidean distance')
plt.ylabel('tRNA score')

plt.savefig(f'./tRNA_vs_5mer_{gene_name}.png')     
plt.close()

print(f"Scatterplot created for {gene_name}")



