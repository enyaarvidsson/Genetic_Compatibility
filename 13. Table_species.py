# Ta ut en gen från euclidean_df_all
# Läs in taxonomy_results_all och filtrera för matchningar (ta bara kolumner vi behöver)
# Beräkna medel-euclidean distance för varje art genen har matchat i (groupby species)
import pandas as pd
import time

start_time = time.time()
gene_name = 'SHV-52'
euclidean_df = pd.read_pickle(f"/storage/enyaa/REVISED/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl") 

#euclidean_df_path = "/storage/enyaa/REVISED/KMER/euclidean_df_all.pkl"
gene_euclidean_df = euclidean_df.T.reset_index()
gene_euclidean_df.columns = ["Bacteria_ID", "Euclidean_distance"]

print(gene_euclidean_df.head())

taxonomy_path = "/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_all.csv"

chunks = []
for chunk in pd.read_csv(taxonomy_path, sep=",", header=None, chunksize=100_000):
    chunk.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    filtered_chunk = chunk[chunk['Gene_name'] == gene_name][['Bacteria_ID', 'Species']]
    chunks.append(filtered_chunk)

gene_taxonomy_df = pd.concat(chunks, ignore_index=True)    
print(gene_taxonomy_df.head())

gene_matching_df = gene_taxonomy_df.merge(gene_euclidean_df, on=['Bacteria_ID'], how="inner")
print(gene_matching_df.head())

mean_species_df = gene_matching_df.groupby('Species').agg(
    Mean_eu = ('Euclidean_distance', 'mean'),
    Counts = ('Bacteria_ID', 'count')
    ).reset_index()

# sortera på storlek - mean
mean_species_df = mean_species_df.sort_values(by="Mean_eu")    
print(mean_species_df.head())

mean_species_df.to_csv(f"/home/enyaa/gene_genome/table_species_{gene_name}.csv")

end_time = time.time()

total_time = (end_time - start_time)/60

print(f"Table created for {gene_name} in {total_time} minutes")

