import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr, wilcoxon

# Original 5mer score
file = "/storage/jolunds/FINAL/5mer_score_table_genes.csv"
table_genes_df = pd.read_csv(file)
table_genes_df['Mean_matches'] = pd.to_numeric(table_genes_df['Mean_matches'], errors='coerce')

# Convert to numeric, force errors (non-numeric values) to NaN
table_genes_df['Mean_matches'] = pd.to_numeric(table_genes_df['Mean_matches'], errors='coerce')
# Drop rows where conversion failed (i.e. where value is NaN)
table_genes_df_filtered = table_genes_df.dropna(subset=['Mean_matches']).copy()

group1_eu = list(table_genes_df_filtered['Mean_matches'])
group2_eu = list(table_genes_df_filtered['Mean'])

stat_eu, p_value_eu = wilcoxon(group1_eu, group2_eu, alternative='less')
print(f"For 5mer score: statistic: {stat_eu}, p-value: {p_value_eu}")

# Euclidean 500 bp
file_500 = "/storage/jolunds/FINAL/length_adjusted_5mer_score_table_genes.csv"
table_genes_500 = pd.read_csv(file_500)

table_genes_500['Mean_matches'] = pd.to_numeric(table_genes_500['Mean_matches'], errors='coerce')

# Convert to numeric, force errors (non-numeric values) to NaN
table_genes_500['Mean_matches'] = pd.to_numeric(table_genes_500['Mean_matches'], errors='coerce')
# Drop rows where conversion failed (i.e. where value is NaN)
table_genes_500_filtered = table_genes_500.dropna(subset=['Mean_matches']).copy()

group1_500 = list(table_genes_500_filtered['Mean_eu_matches'])
group2_500 = list(table_genes_500_filtered['Total_eu_mean'])

stat_500, p_value_500 = wilcoxon(group1_500, group2_500, alternative='less')
print(f"For length-adjusted 5mer score: statistic: {stat_500}, p-value: {p_value_500}")

# tRNA score
file_tRNA = "/storage/jolunds/FINAL/tRNA_table_genes.csv"
table_genes_tRNA = pd.read_csv(file_tRNA)

# Convert to numeric, force errors (non-numeric values) to NaN
table_genes_tRNA['Mean_matches'] = pd.to_numeric(table_genes_tRNA['Mean_matches'], errors='coerce')
# Drop rows where conversion failed (i.e. where value is NaN)
table_genes_tRNA = table_genes_tRNA.dropna(subset=['Mean_matches']).copy()

group1_tRNA = list(table_genes_tRNA['Mean_matches'])
group2_tRNA = list(table_genes_tRNA['Mean'])


stat_tRNA, p_value_tRNA = wilcoxon(group1_tRNA, group2_tRNA, alternative='less')
print(f"For tRNA score: statistic: {stat_tRNA}, p-value: {p_value_tRNA}")

min_eu = table_genes_df['Mean'].min()
min_500 = table_genes_500['Mean'].min()
min = min(min_eu, min_500)

max_eu = table_genes_df['Mean'].max()
max_500 = table_genes_500['Mean'].max()
max = max(max_eu, max_500)

# Scatterplot mean 5mer score vs gene length
plt.figure(figsize=(8, 6))
plt.scatter(table_genes_df_filtered["Mean"], table_genes_df_filtered["Gene_length"], alpha=1, s=10)
plt.xlabel("Mean 5mer score", fontsize=14)
plt.ylabel("Gene length", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(min - 0.002, max + 0.002)

plt.grid(True)
plt.tight_layout()
plt.savefig('./gene_length_5mer_score.png') 
plt.close()

# Scatterplot mean eu vs gene length
plt.figure(figsize=(8, 6))
plt.scatter(table_genes_500["Mean"], table_genes_500["Gene_length"], alpha=1, s=10)
plt.xlabel("Mean length-adjusted 5mer score", fontsize=14)
plt.ylabel("Original gene length", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(min - 0.002, max + 0.002)

plt.grid(True)
plt.tight_layout()
plt.savefig('./gene_length_500bp.png') 
plt.close()

# Scatterplot mean tRNA score vs gene length
plt.figure(figsize=(8, 6))
plt.scatter(table_genes_tRNA["Mean"], table_genes_tRNA["Gene_length"], alpha=1, s=10)
plt.xlabel("Mean tRNA score", fontsize=14)
plt.ylabel("Gene length", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(min - 0.002, max + 0.002)

plt.grid(True)
plt.tight_layout()
plt.savefig('./gene_length_tRNA.png') 
plt.close()

