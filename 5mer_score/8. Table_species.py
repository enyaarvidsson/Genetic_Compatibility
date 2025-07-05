
import pandas as pd
import time
import os

start_time = time.time()

# Choose gene
gene_name = 'SHV-53' 
gene_euclidean_df = pd.read_pickle(f"./FINAL/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl").reset_index(drop=True)

gene_matching_df = gene_euclidean_df[gene_euclidean_df["Match_status"] == "Match"]

mean_species_df = gene_matching_df.groupby('Species').agg(
    Mean_eu = ('Euclidean_distance', 'mean'),
    Counts = ('Bacteria_ID', 'count')
    ).reset_index()

# Sort by euclidean distance
mean_species_df = mean_species_df.sort_values(by="Mean_eu")    

directory = os.path.join('.', 'FINAL', 'TABLE_SPECIES')
os.makedirs(directory, exist_ok=True) 

save_path = f"./FINAL/TABLE_SPECIES/{gene_name}_table_species.csv"
mean_species_df.to_csv(save_path, index=False)

end_time = time.time()
total_time = end_time - start_time
num_species = len(mean_species_df)
print(f"Table created for {gene_name} in {total_time} seconds for {num_species} species")
