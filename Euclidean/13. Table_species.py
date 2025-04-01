# Ta ut en gen från euclidean_df_all
# Läs in taxonomy_results_all och filtrera för matchningar (ta bara kolumner vi behöver)
# Beräkna medel-euclidean distance för varje art genen har matchat i (groupby species)
import pandas as pd
import time

start_time = time.time()

gene_name = 'tet(Q)'
euclidean_df = pd.read_pickle(f"/storage/enyaa/REVISED/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl") 

gene_euclidean_df = euclidean_df.T.reset_index()
gene_euclidean_df.columns = ["Bacteria_ID", "Euclidean_distance"]

taxonomy_path = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/taxonomy_results_{gene_name}.csv"
gene_taxonomy_df = pd.read_csv(taxonomy_path)

gene_matching_df = gene_taxonomy_df.merge(gene_euclidean_df, on=['Bacteria_ID'], how="inner")

mean_species_df = gene_matching_df.groupby('Species').agg(
    Mean_eu = ('Euclidean_distance', 'mean'),
    Counts = ('Bacteria_ID', 'count')
    ).reset_index()

# Sort by euclidean distance
mean_species_df = mean_species_df.sort_values(by="Mean_eu")    

mean_species_df.to_csv(f"/home/jolunds/newtest/2. Table_species/table_species_{gene_name}.csv", index=False)

end_time = time.time()

total_time = end_time - start_time

num_species = len(mean_species_df)
print(f"Table created for {gene_name} in {total_time} seconds for {num_species} species")

