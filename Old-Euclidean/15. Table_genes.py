import pandas as pd
from Bio import SeqIO 
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress

path_bigtable = "/storage/jolunds/REVISED/KMER/big_table.csv"
bigtable_df = pd.read_csv(path_bigtable)

bigtable_df['Gene_name'] = bigtable_df['Gene_name'].replace(r'\?', '/', regex=True)
#print(bigtable_df.tail(20))

phylum_counts = {
    'Pseudomonadota': 1076382,
    'Bacillota': 357869,
    'Campylobacterota': 118184,
    'Actinomycetota': 38108,
    'Bacteroidota': 28352,
    'Spirochaetota': 2905
}
phylum_count_df = pd.DataFrame(list(phylum_counts.items()), columns=['Phylum', 'Count'])

bigtable_df = bigtable_df.merge(phylum_count_df, on='Phylum', how='left')

table_genes_df = bigtable_df.groupby('Gene_name', group_keys=False)[['Mean', 'Count']].apply(
    lambda x: (x['Mean'] * x['Count']).sum() / x['Count'].sum()
).reset_index(name='Total_eu_mean')

filepath = "/storage/enyaa/nucleotide_fasta_protein_homolog_model.fasta"

gene_lengths = {}

for record in SeqIO.parse(filepath, "fasta"):
    header = record.description
    
    positions = header.split("|")[3] 
    start = int(positions.split("-")[0])
    end = int(positions.split("-")[1])
    gene_length = abs(end-start) + 1
    gene_name = header.split("|")[-1].split(" ")[0]
    
    gene_lengths[gene_name] = gene_length
    
gene_lengths_df = pd.DataFrame(list(gene_lengths.items()), columns=['Gene_name', 'Gene_length'])

table_genes_df = table_genes_df.merge(gene_lengths_df, on='Gene_name', how='inner')

# Lägg till
    # minsta värdet
    # minsta matchning + art
    # högsta matchning + art
    # mean för matchningar
    # matchningar i hur många phyla

    

'''
save_path = "/home/enyaa/gene_genome/table_genes_df.csv"
table_genes_df.to_csv(save_path, index=False)

slope, intercept, r_value, p_value, std_err = linregress(table_genes_df['Total_eu_mean'], table_genes_df['Gene_length'])

plt.figure(figsize=(8, 6))
sns.regplot(data=table_genes_df, x='Total_eu_mean', y='Gene_length', lowess=True)
#sns.regplot(x=table_genes_df['Total_eu_mean'], y=table_genes_df['Gene_length'], ci=None)
plt.xlabel("Mean euclidean distance")
plt.ylabel("Gene length")
#plt.title(f"R² = {r_value**2:.3f}")
plt.grid(True)
plt.savefig('/home/enyaa/gene_genome/scatterplot_length_linear.png') 
plt.close()

'''