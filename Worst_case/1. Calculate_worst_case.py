import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

with open ("/storage/enyaa/REVISED/KMER/gene_dist/gene_kmer_distributions.pkl", "rb") as file:
    gene_dictionary = pickle.load(file)
    
gene_name = "SHV-52"

gene_kmer_dist = gene_dictionary.get(gene_name) # Get specific gene

gene_dist_df = pd.DataFrame.from_dict(gene_kmer_dist, orient="index") # Change to dataframe
gene_dist_df.rename(columns={0: gene_name}, inplace=True)
print(gene_dist_df.head())
print(gene_dist_df.sum())

with open ("/storage/enyaa/REVISED/KMER/genome_dist/genome_kmer_distributions_1.pkl", "rb") as file:
    genome_dictionary = pickle.load(file)

#bacteria_id = next(iter(genome_dictionary)) # Get first genome

#genome_kmer_dist = genome_dictionary[bacteria_id]
genome_dist_df = pd.DataFrame.from_dict(genome_dictionary, orient='index').T


#genome_dist_df.rename(columns={0: bacteria_id}, inplace=True)

#print(genome_dist_df.head())

compare_df = gene_dist_df.join(genome_dist_df, how="outer").fillna(0)

#print(compare_df.head())
difference_df = compare_df.sub(compare_df.iloc[:, 0], axis=0).abs()
difference_df = difference_df.drop(columns=[compare_df.columns[0]]) # Drop first column, gene is compared to itself
max_diff = difference_df.max(axis=0)

# Store the result as a DataFrame
worst_case_df = max_diff.to_frame(name="Worst_case").reset_index()
worst_case_df.rename(columns={'index': 'Bacteria_ID'}, inplace=True)

print(worst_case_df.head())

# Add phylum column
path = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_lineage_df = pd.read_csv(path, sep="\t", header=None)
full_lineage_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
phylum_mapping = full_lineage_df[["Bacteria_ID", "Phylum"]] # only the bacteria_id and the respective phylum 

worst_case_df = worst_case_df.merge(phylum_mapping, on='Bacteria_ID', how='left')

print(worst_case_df.head())

# Filter for top phyla
top_phyla = worst_case_df["Phylum"].value_counts().head(6)
worst_case_df = worst_case_df[worst_case_df['Phylum'].isin(top_phyla.index)]
    
# Add match status column
path = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/taxonomy_results_{gene_name}.csv" # this contains matches
    
taxonomy_df = pd.read_csv(path, sep=",") 
   
matching_df = taxonomy_df[['Bacteria_ID']].drop_duplicates()
worst_case_df = worst_case_df.merge(matching_df.assign(Match_status="Match"), on="Bacteria_ID", how="left").fillna("No_match") # here a new column is added to euclidean_gene_df called "Match_status" and it says if there are Match

worst_case_phylum = worst_case_df["Phylum"].nunique()

# REMOVE BACTERIA ID COLUMN
print(len(worst_case_df))
worst_case_df_filtered = worst_case_df[:10_000]
'''
sns.histplot(worst_case_df, x="Worst_case", hue="Match_status", hue_order=["No_match", "Match"],
             multiple="stack", bins=30)
plt.savefig("/home/jolunds/newtest/worst_case.png")
plt.close()
'''

# Plot histogram
g = sns.FacetGrid(worst_case_df_filtered, col="Phylum",  col_order = top_phyla.index, sharey=False, sharex=False,
                  col_wrap = 3) 
g.map_dataframe(sns.histplot, x="Worst_case", hue="Match_status", hue_order=["No_match", "Match"], multiple="stack", bins=10) # orange = match
#g.add_legend(title="Match Status")
g.set_axis_labels("Worst Case", "Number of Bacteria")
g.set_titles(col_template="{col_name}")

plt.savefig("/home/jolunds/newtest/worst_case_try.png")
plt.close()

'''
compare_df['Difference'] = (compare_df.iloc[:,0] - compare_df.iloc[:,1]).abs()
worst_case = compare_df['Difference'].max()
compare_df_sorted = compare_df.sort_values(by="Difference", ascending=False)

print(worst_case)
print(compare_df_sorted.head())
'''


