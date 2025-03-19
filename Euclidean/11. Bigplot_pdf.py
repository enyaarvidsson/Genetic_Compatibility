# Create one big pdf with all the histograms

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np 
import seaborn as sns
import time
import os


start_time = time.time() # Starting time

# GENE_NAMES ----------------------------------------
# Load file with gene_names (sorted)
gene_name_file = "/home/enyaa/gene_genome/gene_names.txt"
gene_names_df = pd.read_csv(gene_name_file, header=None, names=["Gene_name"])
#print(gene_names_df.head())


# TAXONOMY (to get phyla later) ---------------------
# Read full taxonomy file
taxonomy_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_taxonomy_df = pd.read_csv(taxonomy_file, sep="\t", header=None) 
full_taxonomy_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


# PDF -----------------------------------------------
pdfFile = PdfPages("/home/enyaa/gene_genome/bigplot.pdf")


# FOR EACH GENE_NAME --------------------------------
for gene_name in gene_names_df["Gene_name"][:3]: # loops through the gene_names in alphabetical order
#for gene_name in sorted(euclidean_df.index): # loops through the gene_names in alphabetical order
    
    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    # EUCLIDEAN DISTANCE for one gene ---------------
    euclidean_gene_df = pd.read_pickle(f"/storage/enyaa/REVISED/KMER/euclidean_split_genes/euclidean_df_{gene_name}.pkl") 
    #euclidean_gene_df = euclidean_gene_df.melt(var_name="Bacteria_ID", value_name="Euclidean_distance")
    # instead of melt - to make it faster:
    data = euclidean_gene_df.to_numpy()
    bacteria_ids = np.tile(euclidean_gene_df.columns, data.shape[0])
    euclidean_values = data.ravel()
    euclidean_gene_df = pd.DataFrame({"Bacteria_ID": bacteria_ids, "Euclidean_distance": euclidean_values})
        # euclidean_gene_df - has one column Bacteria_ID with the bacteria_ids, and one column with euclidean distance
        # for gene_name
    #print(len(euclidean_gene_df))
    
    # MATCH STATUS (from taxonomy) ------------------
    # Load Taxonomy results
    path = f"/storage/jolunds/REVISED/TAXONOMY/taxonomy_split_genes/taxonomy_results_{gene_name}.csv" # this contains matches
    
    if os.path.exists(path):
        taxonomy_gene_df = pd.read_csv(path, sep=",")

        # Find matching bacteria for this gene
        matching_df = taxonomy_gene_df[['Bacteria_ID']].drop_duplicates() # bacteria_ids that has matched with the gene - UNIQUE
        #matching_df = taxonomy_df[taxonomy_df["Gene_name"] == gene_name][["Bacteria_ID"]] # takes the matching bacteria_id, so we will have the bacteria_id that the gene_name matches with
        #print(len(matching_df))
        #print(matching_df['Bacteria_ID'].nunique())
        
    else: # if taxonomy file doesn't exist - there are no matches for the gene
        matching_df = pd.DataFrame(columns=["Bacteria_ID"]) # create an empty matching_df
        print(f"File not found: {path}")

    # Add match status column using merge
    euclidean_gene_df = euclidean_gene_df.merge(matching_df.assign(Match_status="Match"), on="Bacteria_ID", how="left") # here a new column is added to euclidean_gene_df called "Match_status" and it says if there are Match
    
    euclidean_gene_df["Match_status"] = euclidean_gene_df["Match_status"].fillna("No_match")
        # euclidean_gene_df - has Bacteria_ID, Euclidean_distance and Match_status columns
    #print(len(euclidean_gene_df))

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
    #print(top_phyla)
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
    phylum_counts = phylum_counts.reindex(top_phyla.index) # In same order as top_phyla
    #print(phylum_counts)

    # HISTOGRAM -------------------------------------
    # Create histogram with stacked bars
    g = sns.FacetGrid(downsampled_df, col="Phylum", col_order=phylum_counts.index, sharey=False,
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


