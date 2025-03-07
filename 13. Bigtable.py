# Create a table with average compatiblity for each phylum
import pickle
import pandas as pd
import numpy as np
# Read gene dictionary
with open("/storage/enyaa/REVISED/KMER/gene_dist/gene_kmer_distributions.pkl", "rb") as file: #"rb": read binary
    gene_dictionary = pickle.load(file)

#gene_dictionary_10 = dict(list(gene_dictionary.items())[:10])  
# Read genome dictionary
with open("/storage/enyaa/REVISED/KMER/genome_dist/genome_kmer_distributions_1.pkl", "rb") as file: #"rb": read binary
    genome_dictionary = pickle.load(file)

#genome_dictionary_10k = dict(list(genome_dictionary.items())[:10_000])  

genomes_df = pd.DataFrame.from_dict(genome_dictionary, orient="index").T

# Read taxonomy
taxonomy_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_taxonomy_df = pd.read_csv(taxonomy_file, sep="\t", header=None) 
full_taxonomy_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# Load Taxonomy results
path = "/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_1.csv"
taxonomy_df = pd.read_csv(path, sep=",", header=None)

results = []
for gene_name, kmer_dist in gene_dictionary.items():
    gene_df = pd.DataFrame(list(kmer_dist.items()), columns=["kmer", f"{gene_name}"])
    gene_df.set_index("kmer", inplace=True)
    
    distributions_df = gene_df.join(genomes_df, how="outer")
    distributions_df.fillna(0, inplace=True)
    
    # Calculate euclidean distance
    gene_vector = distributions_df[f"{gene_name}"].values[:, None]
    eu_distances = np.linalg.norm(distributions_df.iloc[:, 1:].values - gene_vector, axis=0)
    
    euclidean_df = pd.DataFrame({
    'Bacteria_ID': distributions_df.columns[1:],  # The genome names (exclude the first column)
    'Euclidean_distance': eu_distances
    })
    
    # Add column for phylum
    euclidean_df = euclidean_df.merge(full_taxonomy_df[["Bacteria_ID", "Phylum"]], on="Bacteria_ID", how="left")
    
    # Filter for the top phyla
    top_phyla = euclidean_df["Phylum"].value_counts().head(4)
    euclidean_df = euclidean_df[euclidean_df["Phylum"].isin(top_phyla.index)]  
    
    # Match / No match, and add to euclidean_df
    filtered_df = taxonomy_df[taxonomy_df.iloc[:,0] == gene_name]   # Takes out the information for specific gene
    filtered_df.columns = ["Gene_ID", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

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
    
    #print(euclidean_df.head())
    
    # Calculate average for each phylum
    average_per_phylum = euclidean_df.groupby('Phylum')['Euclidean_distance'].mean()

    # Count the number of matches per phylum
    matches_per_phylum = euclidean_df[euclidean_df['Match_status'] == 'Match'].groupby('Phylum').size()
    result = {'Gene': gene_name}  # Start with gene name
    #print(matches_per_phylum)
    '''
    # Add average Euclidean distances to the result dictionary
    for phylum, avg_dist in average_per_phylum.items():
        result[f'{phylum}_avg'] = avg_dist
    
    # Add match counts to the result dictionary
    for phylum, match_count in matches_per_phylum.items():
        result[f'{phylum}_matches'] = match_count
    '''
    for phylum, avg_dist in average_per_phylum.items():
        match_count = matches_per_phylum.get(phylum, 0)  # Get match count, defaulting to 0 if no matches
        result[f'{phylum}'] = f"{avg_dist:.3f}, {match_count}"
    
    results.append(result)    


# Convert results list to DataFrame
summary_df = pd.DataFrame(results).fillna(0)  # Fill NaN with 0 for missing values

#print(summary_df.head())
# Save to TSV file
save_path = "/storage/jolunds/big_table.tsv"
summary_df.to_csv(save_path, sep='\t', index=False)    
    