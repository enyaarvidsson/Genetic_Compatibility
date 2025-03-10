# Create a table with average compatiblity for each phylum
import pickle
import pandas as pd
import numpy as np
import time

start_time = time.time()
# Read gene dictionary
with open("/storage/enyaa/REVISED/KMER/gene_dist/gene_kmer_distributions.pkl", "rb") as file: #"rb": read binary
    gene_dictionary = pickle.load(file)

#gene_dictionary_1 = dict(list(gene_dictionary.items())[:10]) #test for 10 genes

# Read genome dictionary
with open("/storage/enyaa/REVISED/KMER/genome_dist/genome_kmer_distributions_1.pkl", "rb") as file: #"rb": read binary
    genome_dictionary = pickle.load(file)

genome_dictionary_10k = dict(list(genome_dictionary.items())[:10_000])  #test for 10k genomes
genomes_df = pd.DataFrame.from_dict(genome_dictionary_10k, orient="index").T # Change tp dataframe

# Read taxonomy
taxonomy_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_taxonomy_df = pd.read_csv(taxonomy_file, sep="\t", header=None) 
full_taxonomy_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# Load Taxonomy results
path = "/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_1.csv" # kör om
taxonomy_df = pd.read_csv(path, sep=",", header=None)
taxonomy_df.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# Count number of matches for each gene in each phylum
phylum_bacteria_counts = taxonomy_df.groupby(["Gene_name", "Phylum"])["Bacteria_ID"].nunique().reset_index()
phylum_bacteria_counts.columns = ["Gene_name", "Phylum", "Unique_Bacteria_Count"] #Change names 

#print(phylum_bacteria_counts)

results = []
for gene_name, kmer_dist in gene_dictionary.items():
    gene_df = pd.DataFrame(kmer_dist.items(), columns=["kmer", f"{gene_name}"]).set_index("kmer")

    # Merge with genome distributions
    distributions_df = gene_df.join(genomes_df, how="outer").fillna(0)
    
    # Calculate euclidean distance
    gene_vector = distributions_df[f"{gene_name}"].to_numpy()[:, np.newaxis]
    eu_distances = np.linalg.norm(distributions_df.iloc[:, 1:].to_numpy() - gene_vector, axis=0)
    
    euclidean_df = pd.DataFrame({
    'Bacteria_ID': distributions_df.columns[1:],  # The genome names (exclude the first column)
    'Euclidean_distance': eu_distances
    })
    
    # Add column for phylum
    euclidean_df = euclidean_df.merge(full_taxonomy_df[["Bacteria_ID", "Phylum"]], on="Bacteria_ID", how="left")
    
    # Filter for the top phyla
    top_phyla = euclidean_df["Phylum"].value_counts().nlargest(4)
    euclidean_df = euclidean_df[euclidean_df["Phylum"].isin(top_phyla.index)]  
    
    # Add gene name
    euclidean_df["Gene_name"] = gene_name
    
    phylum_stats = euclidean_df.merge(phylum_bacteria_counts[["Gene_name", "Phylum", "Unique_Bacteria_Count"]],
                                      on=["Gene_name", "Phylum"], how="left")
    
    

    # Add the mean and standard deviation for each phylum and gene combination
    phylum_stats = phylum_stats.groupby(["Gene_name", "Phylum"]).agg(
        mean=("Euclidean_distance", "mean"),
        std=("Euclidean_distance", "std"),
        Unique_Bacteria_Count=("Unique_Bacteria_Count", "first")  # Get the bacteria count
    ).reset_index()
    
    print(phylum_stats)
    
    results.append(phylum_stats)    

final_df = pd.concat(results, ignore_index=True).fillna(0)
final_df["Unique_Bacteria_Count"] = final_df["Unique_Bacteria_Count"].astype(int)

#print(final_df.head())
# Now you can save it to CSV
final_df.to_csv("/storage/jolunds/big_table_1.csv", index=False)

end_time = time.time()

total_time = (end_time - start_time)/60

print(f"Done creating big table with elapsed time: {total_time} minutes")


'''

summary_df = pd.DataFrame(results).fillna(0)  # Fill NaN with 0 for missing values

#print(summary_df.head())
# Save to TSV file
save_path = "/storage/jolunds/big_table.tsv"
summary_df.to_csv(save_path, sep='\t', index=False)    


    match_column = []
    match = 0
    no_match = 0
    for bacteria in euclidean_df['Bacteria_ID']:
        if bacteria in filtered_df.iloc[:, 1].values:
            match_column.append("Match")
            match += 1
        else:
            match_column.append("No_match")
            no_match += 1

    euclidean_df['Match_status'] = match_column 
    
      # Match / No match, and add to euclidean_df
    #filtered_df = taxonomy_df[taxonomy_df.iloc[:,0] == gene_name]   # Takes out the information for specific gene
    matching_bacteria = taxonomy_df.loc[taxonomy_df["Gene_name"] == gene_name, "Bacteria_ID"].values
    euclidean_df["Match_status"] = np.where(euclidean_df["Bacteria_ID"].isin(matching_bacteria), "Match", "No_match")
    
    # Calculate average and standard eviation for each phylum
    phylum_stats = euclidean_df.groupby("Phylum")["Euclidean_distance"].agg(["mean", "std"])
    
    # Count the number of matches per phylum
    matches_per_phylum = euclidean_df.groupby("Phylum")["Match_status"].apply(lambda x: (x == "Match").sum())
    
    # Store results
    result = {'Gene': gene_name}  
    for phylum in phylum_stats.index:
        avg_dist = phylum_stats.loc[phylum, "mean"]
        std_dist = phylum_stats.loc[phylum, "std"]
        match_count = matches_per_phylum.get(phylum, 0)
        result[f"{phylum}"] = f"{avg_dist:.3f}, {std_dist:.3f}, {match_count}"


    '''