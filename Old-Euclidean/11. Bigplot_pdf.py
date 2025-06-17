# FLYTTAD!

# Create one big pdf with all the histograms

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np 
import seaborn as sns
import time
import os
from tqdm import tqdm
from matplotlib.ticker import MaxNLocator


start_time = time.time() # Starting time

# GENE_NAMES ----------------------------------------
# Load file with gene_names (sorted)
#gene_name_file = "/storage/enyaa/REVISED/gene_names.txt"
#gene_names_df = pd.read_csv(gene_name_file, header=None, names=["Gene_name"])
gene_names_df = pd.DataFrame({"Gene_name": ["tet(Q)"]}) # to pick only one gene


# TAXONOMY (to get phyla later) ---------------------
# Read full taxonomy file
taxonomy_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
dtype_dict = {0: "string", 2: "category"}  
use_cols = [0, 2]  # only load Bacteria_ID and Phylum 
full_taxonomy_df = pd.read_csv(taxonomy_file, sep="\t", usecols=use_cols, dtype=dtype_dict, header=None)
full_taxonomy_df.columns = ["Bacteria_ID", "Phylum"]


# PDF -----------------------------------------------
# to create a pdf - change lines at bottom also
#pdfFile = PdfPages("/home/enyaa/gene_genome/bigplot_500bp.pdf") 


# FOR EACH GENE_NAME --------------------------------
for gene_name in tqdm(gene_names_df["Gene_name"], desc="Processing genes"): # loops through the gene_names in alphabetical order

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    # EUCLIDEAN DISTANCE for one gene ---------------
    euclidean_gene_df = pd.read_pickle(f"/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_{gene_name}.pkl") 
        # euclidean_gene_df - has one column Bacteria_ID with the bacteria_ids, and one column with euclidean distance
    
    # for 500 bp:
    #eu_path = f"/storage/enyaa/REVISED/KMER/FOR_GENE_LENGTH/euclidean_split_genes_500bp/euclidean_df_{gene_name}.pkl"
    #if not os.path.exists(eu_path): # skip genes that are shorter than 500 bp because those files don't exist
    #    continue
    #euclidean_gene_df = pd.read_pickle(eu_path).reset_index(drop=True).T.reset_index() # Switch to long format 
    #euclidean_gene_df.columns = ['Bacteria_ID', 'Euclidean_distance']

    # MATCH STATUS (from taxonomy) ------------------
    # Load Taxonomy results
    path = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/taxonomy_results_{gene_name}.csv" # this contains matches
    
    if os.path.exists(path):
        dtype_dict = {"Bacteria_ID": "string"}
        use_cols = ["Bacteria_ID"]
        taxonomy_gene_df = pd.read_csv(path, sep=",", usecols=use_cols, dtype=dtype_dict)
    
        # Find matching bacteria for this gene
        matching_df = taxonomy_gene_df[['Bacteria_ID']].drop_duplicates() # bacteria_ids that has matched with the gene - UNIQUE
        matches = 1 # there exists matches
        
    else: # if taxonomy file doesn't exist - there are no matches for the gene
        matching_df = pd.DataFrame(columns=["Bacteria_ID"]) # create an empty matching_df
        #print(f"File not found: {path}")
        matches = 0 # no matches

    # Add match status column using merge
    euclidean_gene_df = euclidean_gene_df.merge(matching_df.assign(Match_status="Match"), on="Bacteria_ID", how="left") # here a new column is added to euclidean_gene_df called "Match_status" and it says if there are Match   
    euclidean_gene_df["Match_status"] = euclidean_gene_df["Match_status"].fillna("No_match")
        # euclidean_gene_df - has Bacteria_ID, Euclidean_distance and Match_status columns
    no_match_count = euclidean_gene_df["Match_status"].value_counts().get("Match", 0)
    
    if no_match_count == 0:
        print("No matches for gene:", gene_name)
        matches = 0

    # ADD PHYLUM (from full taxonomy) ---------------
    # Merge with full taxonomy to get Phylum information
    euclidean_gene_df = euclidean_gene_df.merge(full_taxonomy_df[["Bacteria_ID", "Phylum"]], on="Bacteria_ID", how="left")
        # euclidean_gene_df - has Bacteria_ID, Euclidean_distance, Match_status and Phylum columns

    # ONLY PHYLA WITH MOST GENOMES ------------------
    # Select top 6 most frequent phyla
    top_phyla = euclidean_gene_df["Phylum"].value_counts().head(6)
    #print(top_phyla)
    euclidean_top_phyla_df = euclidean_gene_df[euclidean_gene_df["Phylum"].isin(top_phyla.index)]
        # euclidean_top_phyla_df - has only the top 6 phyla - for each gene (columns: Bacteria_ID, Euclidean_distance, Match_status, Phylum)

    # BINS FOR THE HISTOGRAM ------------------------
    nr_bins = 30
    min_value = euclidean_top_phyla_df["Euclidean_distance"].min()
    max_value = euclidean_top_phyla_df["Euclidean_distance"].max() + 0.001 # so all values fall inside the max_value
    bin_edges = np.linspace(min_value, max_value, nr_bins + 1)

    # DOWNSAMPLE NO_MATCH ---------------------------
    # om man vill göra fast för varje phylum och beroende på hur många matches, gör följande
    downsampled_no_matches = [] # will become a list of dataframes

    euclidean_top_phyla_df = euclidean_top_phyla_df.copy()
    euclidean_top_phyla_df["Phylum"] = euclidean_top_phyla_df["Phylum"].astype(str)
    #print(len(euclidean_top_phyla_df))
    #no_match_count = euclidean_top_phyla_df["Match_status"].value_counts().get("No_match", 0)
    #print(no_match_count)

    for phylum, phylum_df in euclidean_top_phyla_df.groupby("Phylum"): # phylum - name of phylum, phylum_df - df that has only rows from that phylum

        match_count = (phylum_df["Match_status"] == "Match").sum()

        # create a df for match and a df for no_match
        matches_phylum_df = phylum_df[phylum_df["Match_status"] == "Match"]
        no_matches_phylum_df = phylum_df[phylum_df["Match_status"] == "No_match"]  

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
        # downsampled_df contains all matches, but downsampled no_matches

    # COUNT NR OF BACTERIA --------------------------
    # in each phylum
    phylum_counts = downsampled_df['Phylum'].value_counts()  
    phylum_counts = phylum_counts.reindex(top_phyla.index) # In same order as top_phyla

    # matches in each phylum
    matches_phylum_counts = downsampled_df[downsampled_df['Match_status'] == 'Match'].groupby('Phylum').size()
    matches_phylum_counts = matches_phylum_counts.reindex(top_phyla.index).fillna(0).astype(int)
    #print(matches_phylum_counts)

    # HISTOGRAM -------------------------------------
    # Create histogram with stacked bars
    g = sns.FacetGrid(downsampled_df, col="Phylum", col_order=phylum_counts.index, sharey=False,
                       col_wrap=3, height=4, aspect=1.2)
    g.map_dataframe(sns.histplot, x="Euclidean_distance", hue="Match_status", hue_order=["No_match", "Match"], 
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
    #if "?" in gene_name:
    #    gene_name = gene_name.replace("?", "/")
    #if matches == 1: # if matches exists
    #    plt.figtext(0.5, 0.95, f"Gene name: {gene_name}", ha="center", fontsize=14)
    #else:
    #    plt.figtext(0.5, 0.95, f"Gene name: {gene_name} - NO MATCHES", ha="center", fontsize=14)     
    
    #pdfFile.savefig(g.figure)
    #plt.close(g.figure)
    plt.tight_layout()
    plt.savefig(f'/home/enyaa/gene_genome/histogram_5mer_{gene_name}.png') # for one gene
    plt.close() # for one gene



# CLOSE PDF -------------------------------------
#pdfFile.close()

print(min_value, "and", max_value)

end_time = time.time()
total_time = (end_time - start_time)/60
print(f"Bigplot_pdf created in: {total_time} minutes")
