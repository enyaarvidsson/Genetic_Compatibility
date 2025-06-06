# Add taxonomic info for the bacteria
# Filter to only keep 
    # 6 largest phyla
    # 10 from each species
# Save in csv file

import pandas as pd
import os


# Get all Bacteria_IDs ---------------
file_paths_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_filepaths.tsv"
df = pd.read_csv(file_paths_file, sep="\t", header=None, names=["Filepath"]) 
    # 1631873 bacteria
df["Filename"] = df["Filepath"].apply(os.path.basename) # gets the filename (not whole filepath) 
df = df.drop(columns=["Filepath"])
df["Bacteria_ID"] = df["Filename"].str.split("_").str[:2].str.join("_") # gets the bacterial ID
df = df.drop(columns=["Filename"])
    # 1631873 bacteria

# Add taxonomic info -----------------
taxonomy_path = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
taxonomy_df = pd.read_csv(taxonomy_path, sep="\t", header=None)
taxonomy_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

all_bacteria_df = df.merge(taxonomy_df, how='left', on='Bacteria_ID')
    # 1631873 bacteria, with taxonomic info

# Keep only 10 of each species -------
species_count = all_bacteria_df["Species"].value_counts()

def filter_species(group):
    if len(group) > 10:
        return group.sample(n=10, random_state=42) 
    return group

filtered_bacteria_df = all_bacteria_df.groupby('Species', group_keys=False).apply(filter_species, include_groups=False).reset_index(drop=True)
filtered_bacteria_df = pd.merge(filtered_bacteria_df, taxonomy_df[['Bacteria_ID', 'Species']], on='Bacteria_ID', how='left')
    # 77182 bacteria, with taxonomic info

# Keep only top 6 phyla --------------
top_phyla = filtered_bacteria_df["Phylum"].value_counts().head(6)
filtered_bacteria_df = filtered_bacteria_df[filtered_bacteria_df["Phylum"].isin(top_phyla.index)]
    # 72690 bacteria, with taxonomic info

# Save to csv file -------------------
path = '/storage/enyaa/FINAL/filtered_bacteria.csv'
filtered_bacteria_df.to_csv(path, index=False)
