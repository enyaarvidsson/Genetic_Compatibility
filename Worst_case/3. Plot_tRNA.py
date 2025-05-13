
# Scatterplot tRNA score vs worst case combined

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler


gene_name = "tet(Q)" 
tRNA_score = "tRNA_score_one_sided"


np.random.seed(42)

if "/" in gene_name:
    gene_name = gene_name.replace("/", "?")


# tRNA SCORE ------------------
filepath_tRNA = f"/storage/jolunds/REVISED/tRNA/tRNA_score_new/tRNA_score_{gene_name}.csv"
tRNA_score_df = pd.read_csv(filepath_tRNA)

# only top phyla
top_phyla = tRNA_score_df["Phylum"].value_counts().head(6)
tRNA_score_df = tRNA_score_df[tRNA_score_df["Phylum"].isin(top_phyla.index)]

# add match status
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
    #print("No matches for gene:", gene_name)
    matches = 0


# DOWNSAMPLE NO-MATCHES -------
tRNA_score_downsampled = [] # will become a list of dataframes

match_count = (tRNA_score_df["Match_status"] == "Match").sum()

# create a df for match and a df for no_match
matches_df = tRNA_score_df[tRNA_score_df["Match_status"] == "Match"]
no_matches_df = tRNA_score_df[tRNA_score_df["Match_status"] == "No_match"]  

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


# WORST CASE ------------------
file_path = f"/storage/jolunds/REVISED/WORST_CASE/worst_case_split/worst_case_{gene_name}.csv"
worst_case_df = pd.read_csv(file_path, sep=",")
worst_case_df = worst_case_df.drop(worst_case_df.columns[0], axis=1)


# MERGE THEM TOGETHER ----------
tRNA_and_worst_df = pd.merge(tRNA_score_downsampled_df, worst_case_df, on='Bacteria_ID', how='inner')

# WEIGHTED COMBINATION OF MAX DIFF AND RELATIVE DIFF -----------
scaler = MinMaxScaler()
tRNA_and_worst_df[["Max_diff_scaled", "Rel_diff_scaled"]] = scaler.fit_transform(
    tRNA_and_worst_df[["Max_difference", "Relative_difference"]]
)

tRNA_and_worst_df["Combined_score"] = (
    0.5 * tRNA_and_worst_df["Max_diff_scaled"] +
    0.5 * tRNA_and_worst_df["Rel_diff_scaled"]
)

# To change name in the legend in the plot
tRNA_and_worst_df['Match_status'] = tRNA_and_worst_df['Match_status'].replace({
    'No_match': 'No match',
    'Match': 'Match'  
})


# SCATTERPLOT -----------
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=tRNA_and_worst_df,
    x='Combined_score',
    y=tRNA_score,   
    hue='Match_status',
    hue_order=["No match", "Match"],
    alpha=1,
    s=20
    #color='darkorange'
)

# Add title
if "?" in gene_name:
    gene_name = gene_name.replace("?", "/")

if tRNA_score == "tRNA_score_one_sided":
    tRNA_score_title = "one-sided"
    tRNA_score_nr = "1"
else:
    tRNA_score_title = "two-sided"
    tRNA_score_nr = "2"
 
#plt.title(f'Matches for {gene_name} ({tRNA_score_title}) - spearman: {correlation:.2f} p={p_value:.2f}')
#plt.title(f'Matches for {gene_name} ({tRNA_score_title})')
plt.xlabel('5mer worst case', fontsize=14)
plt.ylabel('tRNA score', fontsize=14)
plt.tick_params(axis='both', labelsize=12)
plt.legend(fontsize=14, loc="upper center")
plt.tight_layout()

plt.savefig(f'/home/enyaa/gene_genome/tRNA{tRNA_score_nr}_vs_worst_{gene_name}.png')     
plt.close()

print(f"Scatterplot created for {gene_name} {tRNA_score_title}")



