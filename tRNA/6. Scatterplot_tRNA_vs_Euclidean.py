
import matplotlib.pyplot as plt
import pandas as pd
import os
import seaborn as sns
from scipy.stats import spearmanr

# FÖR GAMLA tRNA SCORE !!!

gene_name = "floR" 
tRNA_score = "tRNA_score_one_sided"

if "/" in gene_name:
    gene_name = gene_name.replace("/", "?")


# tRNA SCORE ------------------
filepath_tRNA = f"/storage/jolunds/REVISED/tRNA/tRNA_score/tRNA_score_{gene_name}.csv"
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
    print("No matches for gene:", gene_name)
    matches = 0


# ONLY MATHCES ----------------
tRNA_score_df = tRNA_score_df[tRNA_score_df["Match_status"] == 'Match']


# EUCLIDEAN DISTANCE ----------
filepath_eu = f"/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_{gene_name}.pkl"
euclidean_df = pd.read_pickle(filepath_eu)


# MERGE THEM TOGETHER ----------
tRNA_and_euclidean_df = pd.merge(tRNA_score_df, euclidean_df, on='Bacteria_ID', how='inner')


# SPEARMAN CORRELATION --------
correlation, p_value = spearmanr(tRNA_and_euclidean_df['Euclidean_distance'], tRNA_and_euclidean_df[tRNA_score])

print(f"Spearman correlation: {correlation}")
print(f"P-value: {p_value}")


# SCATTERPLOT -----------
sns.scatterplot(
    data=tRNA_and_euclidean_df,
    x='Euclidean_distance',
    y=tRNA_score,   
    #hue='Match_status',
    s=20,
    color='darkorange'
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
 
plt.title(f'Matches for {gene_name} ({tRNA_score_title}) - spearman: {correlation:.2f} p={p_value:.2f}')
plt.xlabel('Euclidean distance')
plt.ylabel('tRNA score')

plt.savefig(f'/home/enyaa/gene_genome/scatterplot{tRNA_score_nr}_{gene_name}.png')     
plt.close()

print(f"Scatterplot created for {gene_name} {tRNA_score_title}")



