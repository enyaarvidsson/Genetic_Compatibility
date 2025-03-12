# Create one big file with all the plots


# !!!!!!
# Kolla upp hur vi ska stänga figuren
# !!!!!!


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random
import pickle
import pandas as pd
import numpy as np
import seaborn as sns
import time

'''
Innan vi hade euclidean_df

# Read gene dictionary
with open("/storage/enyaa/REVISED/KMER/gene_dist/gene_kmer_distributions.pkl", "rb") as file: #"rb": read binary
    gene_dictionary = pickle.load(file)

gene_dictionary_10 = dict(list(gene_dictionary.items())[:10])

# Read genome dictionary
with open("/storage/enyaa/REVISED/KMER/genome_dist/genome_kmer_distributions_1.pkl", "rb") as file:
    genome_dictionary = pickle.load(file)

# Take 10k random genomes from genome_dictionary:
random.seed(42)
random_genomes = random.sample(list(genome_dictionary.keys()), 10_000)
genome_dictionary_10k = {key: genome_dictionary[key] for key in random_genomes}

# Convert to dataframe 
genomes_df = pd.DataFrame.from_dict(genome_dictionary_10k, orient="index").T
'''


start_time = time.time() # Starting time


# EUCLIDEAN DISTANCE ---------------------------
# Load euclidean_df
distance_df = pd.read_pickle("/storage/enyaa/REVISED/KMER/euclidean_df.pkl") # NOT CREATED YET
    # euclidean_df is a df with gene_name as rows and genome_id as columns, and euclidean distance as values


# ??? ------------------------
# Load Taxonomy results
path = "/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_all.csv"
taxonomy_results_df = pd.read_csv(path, sep=",", header=None)
taxonomy_results_df.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# Read full taxonomy file
taxonomy_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_taxonomy_df = pd.read_csv(taxonomy_file, sep="\t", header=None) 
full_taxonomy_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]



pdfFile = PdfPages("/home/enyaa/gene_genome/bigplot_3.pdf")

for gene_name in distance_df.index[:3]:
    # get the euclidean distances for the gene_name
    euclidean_df = pd.DataFrame({
        'Bacteria_ID': distance_df.columns,
        'Euclidean_distance': distance_df.loc[gene_name]
    })

    # Find matching bacteria for this gene
    matching_df = taxonomy_results_df[taxonomy_results_df["Gene_name"] == gene_name][["Bacteria_ID"]]
    
    # Add match status column using merge
    euclidean_df = euclidean_df.merge(matching_df.assign(Match_status="Match"), on="Bacteria_ID", how="left")
    euclidean_df["Match_status"] = euclidean_df["Match_status"].fillna("No_match")

    # Merge with full taxonomy to get Phylum information
    euclidean_df = euclidean_df.merge(full_taxonomy_df[["Bacteria_ID", "Phylum"]], on="Bacteria_ID", how="left")

    # Select top 6 most frequent phyla
    top_phyla = euclidean_df["Phylum"].value_counts().head(6)
    euclidean_top_phyla_df = euclidean_df[euclidean_df["Phylum"].isin(top_phyla.index)]


    # HISTOGRAM
    nr_bins = 20
    min_value, max_value = euclidean_top_phyla_df["Euclidean_distance"].min(), euclidean_top_phyla_df["Euclidean_distance"].max()
    bin_edges = np.linspace(min_value, max_value, nr_bins + 1)

    # Create histogram with stacked bars
    g = sns.FacetGrid(euclidean_top_phyla_df, col="Phylum", col_order=top_phyla.index, sharey=False,
                       col_wrap=3, height=4, aspect=1.2)
    g.map_dataframe(sns.histplot, x="Euclidean_distance", hue="Match_status", hue_order=["No_match", "Match"], 
                    multiple="stack", bins=bin_edges)
    
    g.set_axis_labels("Euclidean Distance", "Number of Bacteria")

    # **Use `top_phyla` for setting the subplot titles**
    for ax, phylum in zip(g.axes.flat, top_phyla.index):
        ax.set_title(f"{phylum} (n={top_phyla[phylum]})")

    g.set(xlim=(min_value - 0.001, max_value + 0.001))     
    plt.subplots_adjust(top=0.85)

    # Add title
    plt.figtext(0.5, 0.95, f"Gene name: {gene_name}", ha="center", fontsize=14)

    pdfFile.savefig(g.figure)
    plt.close(g.figure)

pdfFile.close()

end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Bigplot_pdf created in: {total_time} minutes")



'''
från innan euclidean_df

for gene_name, kmer_dist in gene_dictionary.items():
   # make dataframe with the gene's kmer distribution
    gene_df = pd.DataFrame(list(kmer_dist.items()), columns=["kmer", f"{gene_name}"])
    gene_df.set_index("kmer", inplace=True)
    
    # merge dataframe with genomes dataframe
    distributions_df = gene_df.join(genomes_df, how="outer")
    distributions_df.fillna(0, inplace=True)
    
    # calculate euclidean distance
    gene_vector = distributions_df[f"{gene_name}"].values[:, None]
    eu_distances = np.linalg.norm(distributions_df.iloc[:, 1:].values - gene_vector, axis=0)
    
    # Store bacteria id with the euclidiean distances 
    euclidean_df = pd.DataFrame({
    'Bacteria_ID': distributions_df.columns[1:],  # The genome names (exclude the first column)
    'Euclidean_distance': eu_distances
    })

    # match, no match
    matching_df = taxonomy_results_df[taxonomy_results_df.iloc[:,0] == gene_name]   # Takes out the information for specific gene
    match_column = []
    for bacteria in euclidean_df['Bacteria_ID']:
        if bacteria in matching_df.iloc[:, 1].values:
            match_column.append("Match")
        else:
            match_column.append("No_match")
    euclidean_df['Match_status'] = match_column # Add column with the match status
    
    # top_phyla
    euclidean_df = euclidean_df.merge(full_taxonomy_df[["Bacteria_ID", "Phylum"]], on="Bacteria_ID", how="left")
    top_phyla = euclidean_df["Phylum"].value_counts().head(6)
    euclidean_top_phyla_df = euclidean_df[euclidean_df["Phylum"].isin(top_phyla.index)]
    
    # Create histogram with stacked bars
    plt.figure(figsize=(8, 5))

    # Make sure the bins are same size for all subplots 
    nr_bins = 20 # Number of bins/bars 
    min_value = euclidean_top_phyla_df["Euclidean_distance"].min()
    max_value = euclidean_top_phyla_df["Euclidean_distance"].max()
    bin_edges = np.linspace(min_value, max_value, nr_bins + 1)

    g = sns.FacetGrid(euclidean_top_phyla_df, col="Phylum",  col_order = top_phyla.index, sharey=False,
                  col_wrap = 3, height=4, aspect=1.2)  
    g.map_dataframe(sns.histplot, x="Euclidean_distance", hue="Match_status", hue_order=["No_match", "Match"], multiple="stack", bins=bin_edges)

    g.set_axis_labels("Euclidean Distance", "Number of Bacteria")

    for ax, phylum in zip(g.axes.flat, top_phyla.index):
        ax.set_title(f"{phylum} (n={top_phyla[phylum]})")

    g.set(xlim=(min_value - 0.001, max_value + 0.001))     
    plt.subplots_adjust(top=0.85)  # Adjust the top margin to make space for the figure title

    #gene_name = gene_name_map.get(gene_id, "Unknown")  # Get gene name, write "Unknown" if not found
    plt.figtext(0.5, 0.95, f"Gene name: {gene_name}", ha="center", fontsize=14)
    
    pdfFile.savefig(g.figure) 
    plt.close(g.figure)
    plt.close()
    
    
    
pdfFile.close()
'''    



