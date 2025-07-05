# Create one big pdf with all the histograms, for the 5mer score
# OR
# create histograms for one specific gene

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np 
import seaborn as sns
import time
import os
from tqdm import tqdm
from matplotlib.ticker import MaxNLocator


start_time = time.time() 

# FOR ALL GENES OR ONE SPECIFIC GENE - change here:
# gene = "all" - for all genes
# gene = "tet(Q)" - for the specific gene tet(Q) - can be changed to any other gene
gene = "all"


# GENE_NAMES ----------------------------------------
if gene == "all":
    # Load file with gene_names (sorted)
    gene_name_file = "./FINAL/gene_names.txt"
    gene_names_df = pd.read_csv(gene_name_file, header=None, names=["Gene_name"])
else:
    gene_names_df = pd.DataFrame({"Gene_name": [gene]}) # to pick only one gene


# PDF -----------------------------------------------
if gene == "all":
    pdfFile = PdfPages("./FINAL/bigplot.pdf") 


# FOR EACH GENE_NAME --------------------------------
for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"):

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    # EUCLIDEAN DISTANCE for one gene ---------------
    euclidean_gene_df = pd.read_pickle(f"./FINAL/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl") 
        # euclidean_gene_df - has Bacteria_ID, Euclidean_distance, Phylum columns and Match_status
    no_match_count = euclidean_gene_df["Match_status"].value_counts().get("Match", 0)
    
    matches = 1 # if there exists matches
    if no_match_count == 0:
        #print("No matches for gene:", gene_name)
        matches = 0 # if no matches exists

    # BINS FOR THE HISTOGRAM ------------------------
    nr_bins = 30
    min_value = euclidean_gene_df["Euclidean_distance"].min()
    max_value = euclidean_gene_df["Euclidean_distance"].max() + 0.001 # so all values fall inside the max_value
    bin_edges = np.linspace(min_value, max_value, nr_bins + 1)

    # DOWNSAMPLE NO_MATCH ---------------------------
    # for each phylum and depending on number of matches
    downsampled_no_matches = [] # will become a list of dataframes

    euclidean_gene_df = euclidean_gene_df.copy()
    euclidean_gene_df["Phylum"] = euclidean_gene_df["Phylum"].astype(str)

    for phylum, phylum_df in euclidean_gene_df.groupby("Phylum"): # phylum - name of phylum, phylum_df - df that has only rows from that phylum

        match_count = (phylum_df["Match_status"] == "Match").sum()

        # create a df for match and a df for no_match
        matches_phylum_df = phylum_df[phylum_df["Match_status"] == "Match"]
        no_matches_phylum_df = phylum_df[phylum_df["Match_status"] == "Non-match"]  

        # How many no_matches to keep
        if matches == 1: # if matches exists
            if match_count == 0:
                keep_size = 10000
            elif match_count < 100:
                keep_size = 2000
            elif match_count < 3000:
                keep_size = 10000
            else:
                keep_size = match_count * 3
        else:
            keep_size = len(phylum_df) 

        # Downsample no_matches
        if keep_size > len(no_matches_phylum_df): # make sure np.random.choice works, so the keep_size is not bigger than the population size
            keep_size = len(no_matches_phylum_df)
        downsampled_no_matches_phylum_df = no_matches_phylum_df.iloc[np.random.choice(len(no_matches_phylum_df), keep_size, replace=False)]

        # Append both "Match" and downsampled "No_match" bacteria
        downsampled_no_matches.append(pd.concat([matches_phylum_df, downsampled_no_matches_phylum_df]))

    # Combine all phyla into the final df
    downsampled_df = pd.concat(downsampled_no_matches, ignore_index=True)
        # downsampled_df contains all matches, but downsampled non_matches

    # COUNT NR OF BACTERIA --------------------------
    # in each phylum
    phylum_counts = downsampled_df['Phylum'].value_counts()
    top_phyla = euclidean_gene_df["Phylum"].value_counts().head(6) # to get the largest phyla first  
    phylum_counts = phylum_counts.reindex(top_phyla.index) # In same order as top_phyla

    # matches in each phylum
    matches_phylum_counts = downsampled_df[downsampled_df['Match_status'] == 'Match'].groupby('Phylum').size()
    matches_phylum_counts = matches_phylum_counts.reindex(top_phyla.index).fillna(0).astype(int)

    # HISTOGRAM -------------------------------------
    # Create histogram with stacked bars
    g = sns.FacetGrid(downsampled_df, col="Phylum", col_order=phylum_counts.index, sharey=False,
                       col_wrap=3, height=4, aspect=1.2)
    g.map_dataframe(sns.histplot, x="Euclidean_distance", hue="Match_status", hue_order=["Non-match", "Match"], 
                    multiple="stack", bins=bin_edges)
    
    g.set_axis_labels("Euclidean distance", "") # Number of bacteria

    for ax, phylum in zip(g.axes.flat, phylum_counts.index):
        ax.set_title(f"{phylum} (n={phylum_counts[phylum]}, m={matches_phylum_counts[phylum]})", fontsize=15)
        ax.set_xlabel("5mer score", fontsize=15)
        #ax.set_ylabel("Number of Bacteria", fontsize=14)
        #ax.xaxis.set_major_locator(MaxNLocator(nbins=8)) # number of ticks on x-axis
        ax.tick_params(axis='both', labelsize=13)
        
    g.set(xlim=(min_value - 0.001, max_value + 0.001))
    
    plt.subplots_adjust(top=0.85)

    # Add title (don't need in report):
    if "?" in gene_name:
        gene_name = gene_name.replace("?", "/")
    if matches == 1: # if matches exists
        plt.figtext(0.5, 0.95, f"Gene name: {gene_name}", ha="center", fontsize=14)
    else:
        plt.figtext(0.5, 0.95, f"Gene name: {gene_name} - NO MATCHES", ha="center", fontsize=14)     
    
    if gene == "all":
        pdfFile.savefig(g.figure)
        plt.close(g.figure)
    else:
        #plt.tight_layout() # in report
        plt.savefig(f'./FINAL/histogram_5mer_{gene_name}.png') # for one gene
        plt.close() # for one gene


if gene == "all":
    pdfFile.close() # close pdf

end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Bigplot_pdf created in: {total_time} minutes")
