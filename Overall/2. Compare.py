import pandas as pd

# tRNA score
tRNA_table_genes = pd.read_csv("./FINAL/tRNA_table_genes.csv")
tRNA_table_genes = tRNA_table_genes[["Gene_name", "Mean", "Num_phyla"]]


# length-adjusted 5mer score
eu500_table_genes = pd.read_csv("./FINAL/length_adjusted_5mer_score_table_genes.csv")
eu500_table_genes = eu500_table_genes[["Gene_name", "Mean", "Num_phyla"]]


# Least compatible 
tRNA_table_genes_least = tRNA_table_genes[tRNA_table_genes["Num_phyla"] >= 3]
tRNA_table_genes_least = tRNA_table_genes_least.sort_values(by="Mean", ascending=False)

eu500_table_genes_least = eu500_table_genes[eu500_table_genes["Num_phyla_match"] >= 3]
eu500_table_genes_least = eu500_table_genes_least.sort_values(by="Mean", ascending=False)

least_tRNA_genes = tRNA_table_genes_least["Gene_name"].head(6)
least_eu500_genes = eu500_table_genes_least["Gene_name"].head(6)

least_compatible = set(least_tRNA_genes).intersection(least_eu500_genes)
print(f"Least compatible genes: {least_compatible}")

# Most compatible
tRNA_table_genes_most = tRNA_table_genes_least.sort_values(by="Mean", ascending=True)

eu500_table_genes_most = eu500_table_genes_least.sort_values(by="Mean", ascending=True)

top10_tRNA_genes = tRNA_table_genes_most["Gene_name"].head(300)
top10_eu500_genes = eu500_table_genes_most["Gene_name"].head(300)

most_compatible = set(top10_tRNA_genes).intersection(top10_eu500_genes)

print(f"Most compatible genes: {most_compatible}")



