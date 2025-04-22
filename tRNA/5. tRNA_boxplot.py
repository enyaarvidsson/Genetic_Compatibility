import pandas as pd
import os
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

# Load one euclidean df file to get the filtered bacteria ids (77k)
file = "/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_tet(Q).pkl"
filtered_df = pd.read_pickle(file)
bacteria_ids = filtered_df[["Bacteria_ID"]]


# Load full lineage and add phylum and species info 
path = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_lineage_df = pd.read_csv(path, sep="\t", header=None)
full_lineage_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
phylum_mapping = full_lineage_df[["Bacteria_ID", "Phylum"]] # only the bacteria_id and the respective phylum 

bacteria_ids_phylum = bacteria_ids.merge(phylum_mapping, on="Bacteria_ID", how="left")

from collections import Counter

all_counts = []

for bacteria_id in bacteria_ids["Bacteria_ID"]:
    tRNA_file = f"/storage/jolunds/REVISED/tRNA/tRNA_results/{bacteria_id}_trnascan.txt"
    if not os.path.exists(tRNA_file):
        continue
    
    tRNA_df = pd.read_csv(
        tRNA_file, sep="\t", comment="#", skiprows=3, header=None,
        names=["Sequence_name", "tRNA_nr", "Begin", "End", "tRNA_type",
               "Anticodon", "Intron_begin", "Intron_end", "Score", "Comment"]
    )
    
    counts = Counter(tRNA_df["Anticodon"])
    counts["Bacteria_ID"] = bacteria_id
    all_counts.append(counts)

# Create DataFrame with anticodons as columns
anticodon_counts_df = pd.DataFrame(all_counts).fillna(0)

# Merge with bacteria_ids_phylum
anticodon_counts_df = bacteria_ids_phylum.merge(anticodon_counts_df, on="Bacteria_ID", how="right").groupby("Phylum")

top_phyla = ["Actinomycetota", "Bacillota", "Bacteroidota", "Campylobacterota", "Cyanobacteriota", "Pseudomonadota"]

anticodon_counts_df = anticodon_counts_df[anticodon_counts_df["Phylum"].isin(top_phyla)]

print(anticodon_counts_df.head())

'''# Boxplot
phylum = "Actinomycetota"

plt.figure(figsize=(14, 6))
plt.boxplot(phylum_anticodon_means.loc[phylum].values, showfliers=False)
plt.xticks(rotation=90)

plt.savefig(f"/home/jolunds/newtest/boxplot_{phylum}.png")


print(len(phylum_anticodon_means.loc[phylum].values))
print(len(anticodon_columns))'''
