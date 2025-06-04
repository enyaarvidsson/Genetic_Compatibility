
# Scatterplot Euclidean distance vs worst case combined

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler


gene_name = "tet(Q)" 


np.random.seed(42)

if "/" in gene_name:
    gene_name = gene_name.replace("/", "?")


# EUCLIDEAN DISTANCE ------------------
filepath_eu = f"/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_{gene_name}.pkl"
euclidean_df = pd.read_pickle(filepath_eu)
    # 77182 rows

# add match status
taxonomy_file = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/taxonomy_results_{gene_name}.csv"
if os.path.exists(taxonomy_file):
    taxonomy_df = pd.read_csv(taxonomy_file)
    taxonomy_df.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    taxonomy_df = taxonomy_df.drop_duplicates() # Remove duplicates (when we have multiple matches in one genome)
    matching_df = taxonomy_df[['Bacteria_ID']]
    euclidean_df = euclidean_df.merge(matching_df.assign(Match_status="Match"), on="Bacteria_ID", how="left").fillna("No_match") # here a new column is added to euclidean_gene_df called "Match_status" and it says if there are Match
    matches = 1
else: # If taxonomy file do not exist
    euclidean_df["Match_status"] = "No_match"
    matches = 0

no_match_count = euclidean_df["Match_status"].value_counts().get("Match", 0)    
if no_match_count == 0:
    #print("No matches for gene:", gene_name)
    matches = 0

# add phyla - from full taxonomy file
taxonomy_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
dtype_dict = {0: "string", 2: "category"}  
use_cols = [0, 2]  # only load Bacteria_ID and Phylum 
full_taxonomy_df = pd.read_csv(taxonomy_file, sep="\t", usecols=use_cols, dtype=dtype_dict, header=None)
full_taxonomy_df.columns = ["Bacteria_ID", "Phylum"]

euclidean_df = euclidean_df.merge(full_taxonomy_df, on="Bacteria_ID", how="left")

# only top phyla
top_phyla = euclidean_df["Phylum"].value_counts().head(6)
euclidean_df = euclidean_df[euclidean_df["Phylum"].isin(top_phyla.index)]
    # 72690 rows


# DOWNSAMPLE NO-MATCHES -------
euclidean_downsampled = [] # will become a list of dataframes

match_count = (euclidean_df["Match_status"] == "Match").sum()

# create a df for match and a df for no_match
matches_df = euclidean_df[euclidean_df["Match_status"] == "Match"]
no_matches_df = euclidean_df[euclidean_df["Match_status"] == "No_match"]  

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
euclidean_downsampled.append(pd.concat([matches_df, downsampled_no_matches_df]))

euclidean_downsampled_df = pd.concat(euclidean_downsampled, ignore_index=True)


# WORST CASE ------------------
file_path = f"/storage/jolunds/REVISED/WORST_CASE/worst_case_split/worst_case_{gene_name}.csv"
worst_case_df = pd.read_csv(file_path, sep=",")
worst_case_df = worst_case_df.drop(worst_case_df.columns[0], axis=1)


# MERGE THEM TOGETHER ----------
euclidean_and_worst_df = pd.merge(euclidean_downsampled_df, worst_case_df, on='Bacteria_ID', how='inner')
#print(euclidean_and_worst_df)

# WEIGHTED COMBINATION OF MAX DIFF AND RELATIVE DIFF -----------
scaler = MinMaxScaler()
euclidean_and_worst_df[["Max_diff_scaled", "Rel_diff_scaled"]] = scaler.fit_transform(
    euclidean_and_worst_df[["Max_difference", "Relative_difference"]]
)
#print(tRNA_and_worst_df)

euclidean_and_worst_df["Combined_score"] = (
    0.5 * euclidean_and_worst_df["Max_diff_scaled"] +
    0.5 * euclidean_and_worst_df["Rel_diff_scaled"]
)

# To change name in the legend in the plot
euclidean_and_worst_df['Match_status'] = euclidean_and_worst_df['Match_status'].replace({
    'No_match': 'Non-match',
    'Match': 'Match'  
})

# SCATTERPLOT -----------
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=euclidean_and_worst_df,
    x='Combined_score',
    y='Euclidean_distance',   
    hue='Match_status',
    hue_order=["Non-match", "Match"],
    alpha=1,
    s=20
    #color='darkorange'
)

# Add title
if "?" in gene_name:
    gene_name = gene_name.replace("?", "/")
 
#plt.title(f'Matches for {gene_name} ({tRNA_score_title}) - spearman: {correlation:.2f} p={p_value:.2f}')
#plt.title(f'{gene_name}')
plt.xlabel('5mer worst case', fontsize=14)
plt.ylabel('5mer score', fontsize=14)
plt.tick_params(axis='both', labelsize=12)
plt.legend(fontsize=14, loc="upper center")
plt.tight_layout()

plt.savefig(f'/home/enyaa/gene_genome/Euclidean_vs_worst_{gene_name}.png')     
plt.close()

print(f"Scatterplot created for {gene_name}")



