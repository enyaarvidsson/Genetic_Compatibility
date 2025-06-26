import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import os

# Load reference data from 5mer score
file = f"/storage/jolunds/FINAL/REFERENCES/5mer_score_reference_df.csv"
reference_5mer_df = pd.read_csv(file)

reference_genes = reference_5mer_df["Gene_name", "Bacteria_ID", "Reference"]

# Loop through length-adjusted 5mer score
reference_list = []

for gene_name in reference_genes["Gene_name"].unique():
    gene_file = f"scores_folder/{gene_name}_scores.csv"  
    
    if os.path.exists(gene_file):
        gene_df = pd.read_csv(gene_file)
        
        # Subset the reference df for this gene
        ref_subset = reference_genes[reference_genes["Gene_name"] == gene_name]

        # Merge on Gene_name and Bacteria_ID
        merged = pd.merge(ref_subset, gene_df, on=["Gene_name", "Bacteria_ID"], how="inner")

        reference_list.append(merged)

# Concat
reference_df = pd.concat(reference_list, ignore_index=True)

compatible_df = reference_df[reference_df["Reference"] == "Compatible"]
incompatible_df = reference_df[reference_df["Reference"] == "Incompatible"]

# Split van & beta
van_incompatible_df = incompatible_df[incompatible_df["Phylum"] == "Pseudomonadota"]
beta_incompatible_df = incompatible_df[incompatible_df["Phylum"] == "Bacillota"]

group1 = list(compatible_df['tRNA_score'])
group2_van = list(van_incompatible_df['tRNA_score'])
group2_beta = list(beta_incompatible_df['tRNA_score'])

stat_van, p_value_van = mannwhitneyu(group1, group2_van, alternative='less')  # Test if group1 has smaller values
stat_beta, p_value_beta = mannwhitneyu(group1, group2_beta, alternative='less') 

reference_van = pd.concat([compatible_df, van_incompatible_df], ignore_index=True)
reference_beta = pd.concat([compatible_df, beta_incompatible_df], ignore_index=True)

# Plot reference genes - tRNA score
plt.figure(figsize=(8, 6))
ax =sns.histplot(data=reference_van, x='tRNA_score', hue='Reference', hue_order=["Compatible", "Incompatible"], multiple='stack', bins=30, palette={'Compatible': 'mediumseagreen', 'Incompatible': 'palevioletred'})

plt.setp(ax.get_legend().get_texts(), fontsize='16')
ax.legend_.set_title(None)
plt.xlabel("tRNA score", fontsize=16)
plt.ylabel("", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig("/storage/jolunds/FINAL/REFERENCES/tRNA_score_referneces_van.png")
plt.close()

print(f"p-value for tRNA score with van is {p_value_van}")

plt.figure(figsize=(8, 6))
ax = sns.histplot(data=reference_beta, x='tRNA_score', hue='Reference', hue_order=["Compatible", "Incompatible"], multiple='stack', bins=30, palette={'Compatible': 'mediumseagreen', 'Incompatible': 'palevioletred'})

plt.setp(ax.get_legend().get_texts(), fontsize='16')
ax.legend_.set_title(None)
plt.xlabel("tRNA score", fontsize=16)
plt.ylabel("", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig("/storage/jolunds/FINAL/REFERENCES/tRNA_score_referneces_beta.png")
plt.close()

print(f"p-value for tRNA score with beta-lactamases is {p_value_beta}")