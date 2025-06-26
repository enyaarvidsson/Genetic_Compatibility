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

group1 = list(compatible_df['Euclidean_distance'])
group2 = list(incompatible_df['Euclidean_distance'])
stat, p_value = mannwhitneyu(group1, group2, alternative='less')  # Test if group1 has smaller values

# Plot reference genes - length adjusted 5mer score
plt.figure(figsize=(8, 6))
g = sns.histplot(data=reference_df, x='Euclidean_distance', hue='Reference', hue_order=["Compatible", "Incompatible"], multiple='stack', bins=30, palette={'Compatible': 'mediumseagreen', 'Incompatible': 'palevioletred'})

plt.setp(g.get_legend().get_texts(), fontsize=16)
g.legend_.set_title(None)
plt.xlabel("Length-adjusted 5mer score", fontsize=16)
plt.ylabel("", fontsize=16)
plt.tight_layout()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig("/storage/jolunds/FINAL/REFERENCES/length_adjusted_5mer_score_references.png")
plt.close()

print(f"p-value for length-adjusted 5mer score: {p_value}")