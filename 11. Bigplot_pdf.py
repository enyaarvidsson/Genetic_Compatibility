# Create one big pdf with all the histograms

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random
import pickle
import pandas as pd
import numpy as np
import seaborn as sns
import time
import os

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

# GENE_NAMES ----------------------------------------
# Load file with gene_names (sorted)
gene_name_file = "/home/enyaa/gene_genome/gene_names.txt"
gene_names_df = pd.read_csv(gene_name_file, header=None, names=["Gene_name"])
#print(gene_names_df.head())

'''
# EUCLIDEAN DISTANCE --------------------------------
# Load euclidean_df
euclidean_df = pd.read_pickle("/storage/enyaa/REVISED/KMER/euclidean_df.pkl") 
    # euclidean_df is a df with gene_name as rows and genome_id as columns, and euclidean distance as values


# TAXONOMY (to find matches later) ------------------
# Load Taxonomy results
path = "/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_1.csv" # this contains matches
taxonomy_df = pd.read_csv(path, sep=",", header=None)
taxonomy_df.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
'''

# TAXONOMY (to get phyla later) ---------------------
# Read full taxonomy file
taxonomy_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_taxonomy_df = pd.read_csv(taxonomy_file, sep="\t", header=None) 
full_taxonomy_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

'''
# SAME X-ASIS FOR ALL THE PLOTS ---------------------
global_min = euclidean_df.min().min() # min x-value
global_max = euclidean_df.max().max() # max x-value
'''

# PDF -----------------------------------------------
pdfFile = PdfPages("/home/enyaa/gene_genome/bigplot_3_all.pdf")


# FOR EACH GENE_NAME --------------------------------
for gene_name in gene_names_df["Gene_name"]: # loops through the gene_names in alphabetical order
#for gene_name in sorted(euclidean_df.index): # loops through the gene_names in alphabetical order
    
    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    # EUCLIDEAN DISTANCE for one gene ---------------
    euclidean_gene_df = pd.read_pickle(f"/storage/enyaa/REVISED/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl") 
    euclidean_gene_df = euclidean_gene_df.melt(var_name="Bacteria_ID", value_name="Euclidean_distance")
        # euclidean_gene_df - has one column Bacteria_ID with the bacteria_ids, and one column with euclidean distance
        # for gene_name
    '''
    euclidean_gene_df = pd.DataFrame({
        'Bacteria_ID': euclidean_df.columns,
        'Euclidean_distance': euclidean_df.loc[gene_name]
    })  # euclidean_gene_df - has one column Bacteria_ID with the bacteria_ids, and one column with euclidean distance (also has bacteria_id as row index)
        # for gene_name 
    '''
    
    # MATCH STATUS (from taxonomy) ------------------
    # Load Taxonomy results
    path = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/taxonomy_results_{gene_name}.csv" # this contains matches
    
    # Check if file exists before reading
    if not os.path.exists(path):
        print(f"File not found: {path}, skipping...")
        continue  # Skip to the next gene_name

    taxonomy_gene_df = pd.read_csv(path, sep=",")
    #print(taxonomy_gene_df)

    # Find matching bacteria for this gene
    matching_df = taxonomy_gene_df[['Bacteria_ID']] # bacteria_ids that has matched with the gene
    #matching_df = taxonomy_df[taxonomy_df["Gene_name"] == gene_name][["Bacteria_ID"]] # takes the matching bacteria_id, so we will have the bacteria_id that the gene_name matches with

    # Add match status column using merge
    euclidean_gene_df = euclidean_gene_df.merge(matching_df.assign(Match_status="Match"), on="Bacteria_ID", how="left") # here a new column is added to euclidean_gene_df called "Match_status" and it says if there are Match
    euclidean_gene_df["Match_status"] = euclidean_gene_df["Match_status"].fillna("No_match")
        # euclidean_gene_df - has Bacteria_ID, Euclidean_distance and Match_status columns

    # ADD PHYLUM (from full taxonomy) ---------------
    # Merge with full taxonomy to get Phylum information
    euclidean_gene_df = euclidean_gene_df.merge(full_taxonomy_df[["Bacteria_ID", "Phylum"]], on="Bacteria_ID", how="left")
        # euclidean_gene_df - has Bacteria_ID, Euclidean_distance, Match_status and Phylum columns

    # ONLY PHYLA WITH MOST GENOMES ------------------
    # Select top 6 most frequent phyla
    top_phyla = euclidean_gene_df["Phylum"].value_counts().head(6)
    euclidean_top_phyla_df = euclidean_gene_df[euclidean_gene_df["Phylum"].isin(top_phyla.index)]
        # euclidean_top_phyla_df - has only the top 6 phyla - for each gene
    #print(euclidean_top_phyla_df)

    # BINS FOR THE HISTOGRAM ------------------------
    nr_bins = 30
    #bin_edges = np.linspace(global_min, global_max, nr_bins + 1)
    min_value = euclidean_top_phyla_df["Euclidean_distance"].min()
    max_value = euclidean_top_phyla_df["Euclidean_distance"].max() + 0.001 # so all values fall inside the max_value
    #min_value, max_value = euclidean_top_phyla_df["Euclidean_distance"].min(), euclidean_top_phyla_df["Euclidean_distance"].max()
    bin_edges = np.linspace(min_value, max_value, nr_bins + 1)

    # DOWNSAMPLE NO_MATCH ---------------------------
    downsample_factor = 0.2 # keep 20% of the no_match bacteria, from each bin

    downsampled_no_matches = [] # will become a list of dataframes

    for phylum, phylum_df in euclidean_top_phyla_df.groupby("Phylum"): # phylum - name of phylum, phylum_df - df that has only rows from that phylum
    
        # create a df for match and a df for no_match
        matches_phylum_df = phylum_df[phylum_df["Match_status"] == "Match"]
        no_matches_phylum_df = phylum_df[phylum_df["Match_status"] == "No_match"]  

        # each bacteria is assigned to a bin
        bin_nr = np.digitize(no_matches_phylum_df["Euclidean_distance"], bin_edges) - 1  # subtracting 1 to get index starting from 0
        no_matches_phylum_df = no_matches_phylum_df.copy()
        no_matches_phylum_df.loc[:, "bin_nr"] = bin_nr
        #no_matches_phylum_df["bin_nr"] = bin_nr

        # group by bins and downsample
        downsampled_no_matches_phylum_df = (
            no_matches_phylum_df.groupby(bin_nr, group_keys=False) # group_keys=False - group labels not added to the output
            .apply(lambda x: x.sample(frac=downsample_factor) if len(x) > 1 else x) 
        ) # lambda is a funtion, x - one bin from groupby, sample - randomly selects rows, if it is only 1 row/bacteria, it keeps it

        # Append both "Match" and downsampled "No_match" bacteria
        downsampled_no_matches.append(pd.concat([matches_phylum_df, downsampled_no_matches_phylum_df]))

    # Combine all phyla into the final df
    downsampled_df = pd.concat(downsampled_no_matches, ignore_index=True)
        # downsampled_df contains all matches, but downsampled no_matches
    #print(downsampled_df.head(10))

    # COUNT NR OF BACTERIA --------------------------
    # in each phylum
    phylum_counts = downsampled_df['Phylum'].value_counts() 
    #print(phylum_counts)

    # HISTOGRAM -------------------------------------
    # Create histogram with stacked bars
    g = sns.FacetGrid(downsampled_df, col="Phylum", col_order=top_phyla.index, sharey=False,
                       col_wrap=3, height=4, aspect=1.2)
    g.map_dataframe(sns.histplot, x="Euclidean_distance", hue="Match_status", hue_order=["No_match", "Match"], 
                    multiple="stack", bins=bin_edges)
    
    g.set_axis_labels("Euclidean Distance", "Number of Bacteria")

    # **Use `top_phyla` for setting the subplot titles**
    for ax, phylum in zip(g.axes.flat, phylum_counts.index):
        ax.set_title(f"{phylum} (n={phylum_counts[phylum]})")
        
    g.set(xlim=(min_value - 0.001, max_value + 0.001))
    
    plt.subplots_adjust(top=0.85)

    # Add title
    if "?" in gene_name:
        gene_name = gene_name.replace("?", "/")
    plt.figtext(0.5, 0.95, f"Gene name: {gene_name}", ha="center", fontsize=14)

    pdfFile.savefig(g.figure)
    plt.close(g.figure)


# CLOSE PDF -------------------------------------
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
    matching_df = taxonomy_df[taxonomy_df.iloc[:,0] == gene_name]   # Takes out the information for specific gene
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



