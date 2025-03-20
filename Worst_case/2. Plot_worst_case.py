import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

gene_name = "ErmB"

file_path = f"/storage/jolunds/REVISED/WORST_CASE/worst_case_split/worst_case_{gene_name}.csv"
worst_case_df = pd.read_csv(file_path, sep=",").T.reset_index()
worst_case_df.rename(columns={'index': 'Bacteria_ID'})

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

# Ta bort hälften av no_match (no_match: 94302)
no_match_count = (worst_case_df['Match_status'] == 'no_match').sum()
no_match_indices = worst_case_df[worst_case_df['Match_status'] == 'No_match'].index

num_to_remove = no_match_count**0.8
random_indices = np.random.choice(no_match_indices, num_to_remove, replace=False)

worst_case_df_filtered= worst_case_df.drop(random_indices).reset_index(drop=True)

print(len(worst_case_df_filtered))
#worst_case_df_filtered = worst_case_df[:20_000]

sns.histplot(worst_case_df, x="Worst_case", hue="Match_status", hue_order=["No_match", "Match"],
             multiple="stack", bins=30)
plt.savefig("/home/jolunds/newtest/worst_case.png")
plt.close()

x_min = worst_case_df_filtered["Worst_case"].min()
x_max = worst_case_df_filtered["Worst_case"].max()
num_bins = 15  # Set the desired number of bins
bin_edges = np.linspace(x_min, x_max, num_bins + 1) 
# Plot histogram
g = sns.FacetGrid(worst_case_df_filtered, col="Phylum",  col_order = top_phyla.index, sharey=False,
                  col_wrap = 3) 
g.map_dataframe(sns.histplot, x="Worst_case", hue="Match_status", hue_order=["No_match", "Match"], multiple="stack", bins=bin_edges) # orange = match
#g.add_legend(title="Match Status")
g.set_axis_labels("Worst Case", "Number of Bacteria")
g.set_titles(col_template="{col_name}")

plt.savefig(f"/home/jolunds/newtest/worst_case{gene_name}.png")
plt.close()



