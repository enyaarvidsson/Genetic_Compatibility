# Create a table with average compatiblity for each phylum
import pandas as pd
import time
import matplotlib.pyplot as plt
#import seaborn as sns

############ ALLA TAXONOMY RESULTS EXISTERAR INTE #########

start_time = time.time()

# Load gene names
with open("/storage/jolunds/REVISED/gene_names.txt", "r") as f:
    all_genes = [line.strip() for line in f]
    #first_gene = f.readline().strip()

# Load full lineage 
path = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_lineage_df = pd.read_csv(path, sep="\t", header=None)
full_lineage_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

phylum_mapping = full_lineage_df[["Bacteria_ID", "Phylum"]] # only the bacteria_id and the respective phylum 

# Find for top phylum
top_phyla = phylum_mapping.groupby('Phylum').size().sort_values(ascending=False).head(6)
##### FÅR EJ SAMMA TOP PHYLA SOM I BIGPLOT 
    # Chloroflexota  ist för Spirochaetota

# Loop through genes in all_genes 
big_table_list = []   
for gene_name in sorted(all_genes): 
    
    if "/" in gene_name: # Just for look-up
        gene_name = gene_name.replace("/", "?")
    
    # Load euclidean df
    gene_euclidean_df = pd.read_pickle(f"/storage/enyaa/REVISED/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl")
    
    # Switch to long format 
    gene_euclidean_df = gene_euclidean_df.melt(ignore_index=False, var_name="Bacteria_ID", value_name="Euclidean_distance").reset_index()
    
    # Add phylum column
    phylum_euclidean_df = gene_euclidean_df.merge(phylum_mapping, on=['Bacteria_ID'], how='inner')
    phylum_euclidean_df.rename(columns={'index': 'Gene_name'})
    
    # Filter for top phyla
    phylum_euclidean_df = phylum_euclidean_df[phylum_euclidean_df['Phylum'].isin(top_phyla.index.tolist())]
    
    # Calculate mean and standard deviation for each phylum and gene
    phylum_stats = phylum_euclidean_df.groupby(["Phylum"])["Euclidean_distance"].agg(['mean', 'std']).fillna(0)
    phylum_stats_reset = phylum_stats.reset_index()
        # has the genes as rows, and shows mean and std for each phylum
    
    # Load taxonomy results and count matches
    taxonomy_df = pd.read_csv(f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/taxonomy_results_{gene_name}.csv")
    taxonomy_df.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    phylum_counts = taxonomy_df.groupby('Phylum').size().reset_index(name="Num_matches")
    
    # Merge stats and num matches
    gene_table_df = pd.merge(phylum_stats_reset, phylum_counts, on=["Phylum"], how="left").fillna(0) # Add num matches 
    
    # Changes for the eye
    gene_table_df.insert(0, 'Gene_name', gene_name)
    gene_table_df.rename(columns={'mean': 'Mean', 'std': 'Std'}, inplace=True)
    gene_table_df['Num_matches'] = gene_table_df['Num_matches'].astype(int)
    gene_table_df['Phylum'] = pd.Categorical(gene_table_df['Phylum'], categories=top_phyla.index, ordered=True)
    gene_table_df = gene_table_df.sort_values('Phylum')
    
    # Append results
    big_table_list.append(gene_table_df)

# Concat
big_table_df = pd.concat(big_table_list).reset_index(drop=True)

# Save
save_path = "/storage/enyaa/REVISED/KMER/big_table.csv"
big_table_df.to_csv(save_path, index=False)

'''
# Create scatterplot
# Använda facet grid för att göra en figur med en plot för varje phylum
g = sns.FacetGrid(big_table_df, col="Phylum", col_order = top_phyla.index, col_wrap = 3)
g.map_dataframe(sns.scatterplot, x="Num_matches", y="Mean")
g.set_axis_labels("Mean euclidean distance", "Number of bacteria")

# Save plot
plt.savefig("/home/jolunds/newtest/scatterplot_test.png")
plt.close(g.figure)
end_time = time.time()
total_time = (end_time - start_time)/60

print(f"Done creating big table and scatterplot with elapsed time: {total_time} minutes")
'''


# GAMMAL VERSION:
'''
euclidean_df = pd.read_pickle("/storage/enyaa/REVISED/KMER/euclidean_df.pkl") # NOT CREATED YET

# Switch to long format 
euclidean_df= euclidean_df.melt(ignore_index=False, var_name="Bacteria_ID", value_name="Euclidean_Distance").reset_index()

# Add phylum column
path = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_lineage_df = pd.read_csv(path, sep="\t", header=None)
full_lineage_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

phylum_mapping = full_lineage_df[["Bacteria_ID", "Phylum"]] # only the bacteria_id and the respective phylum 
phylum_euclidean_df = euclidean_df.merge(phylum_mapping, left_on="Bacteria_ID", right_on="Bacteria_ID")
    # phylum_euclidean_df - genes, genomes, phylum

# Calculate mean and standard deviation for each phylum and gene
phylum_stats = phylum_euclidean_df.groupby(["index", "Phylum"])["Euclidean_Distance"].agg(['mean', 'std']).fillna(0)
phylum_stats_reset = phylum_stats.reset_index()
phylum_stats_reset.rename(columns={'index': 'Gene_name'}, inplace=True)
    # has the genes as rows, and shows mean and std for each phylum

print(phylum_stats_reset.head(10))
# Count number of matches for each gene in each phylum
taxonomy_df = pd.read_csv("/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_all.csv", header=None) # create a pandas dataframe
taxonomy_df.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"] # add column names

phylum_bacteria_counts = taxonomy_df.groupby(["Gene_name", "Phylum"])["Bacteria_ID"].nunique().reset_index()
phylum_bacteria_counts.columns = ["Gene_name", "Phylum", "Num_matches"] 

print(phylum_bacteria_counts.head(10))
# !!!!!!!! HÄR BLIR DET FEL !!!!!!!!!!!!! 
#phylum_stats_reset = phylum_stats.reset_index() 
#big_table_df = pd.merge(phylum_stats_reset, phylum_bacteria_counts, on=["Gene_name", "Phylum"], how="left") # Add num matches 

# BYT TILL how='left' NÄR VI HAR euclidean_df-all
big_table_df = phylum_stats_reset.merge(phylum_bacteria_counts, on=["Gene_name", "Phylum"], how="inner")

print(big_table_df.head(10))
# Filter for the top phyla
phylum_total_counts = phylum_bacteria_counts.groupby("Phylum")["Num_matches"].sum().reset_index()
print(phylum_total_counts)
top_phyla = phylum_total_counts.sort_values(by="Num_matches", ascending=False).head(6)["Phylum"]

big_table_df = big_table_df[big_table_df["Phylum"].isin(top_phyla)]
print(big_table_df.head(10))
'''
