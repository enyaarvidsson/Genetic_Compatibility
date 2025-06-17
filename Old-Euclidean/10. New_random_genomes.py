# TROR INTE DENNA BEHÖVS?

# New random genomes

# Takes one gene, and many genomes
# It calculates the eucledian distance between the gene and all of the genomes
# 


# Gene 1: AE006468.2, mobile, 409
# Gene 2: Y10279.1,mobile,514
import pickle
import random 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# GENE -------------------------------------------
# Read gene dictionary
with open("/storage/enyaa/REVISED/KMER/gene_dist/gene_kmer_distributions.pkl", "rb") as file: # "rb": read binary
    gene_dictionary = pickle.load(file) # for each gene name, we have the kmers (ex AAAAC) and the corresponding values (ex 0.0043) 

# Find the desired gene
gene_name = "tet(Q)" # Gene name for a gene that we picked
gene_distribution = gene_dictionary.get(gene_name) # Get the desired gene dist from the dictionary


# GENOMES -----------------------------------------
# Read genome dictionary
with open("/storage/enyaa/REVISED/KMER/genome_dist/genome_kmer_distributions_1.pkl", "rb") as file:
    genome_dictionary = pickle.load(file)

# Find random genomes
random.seed(50) # get the same genomes each time
random_bacteria_ids = random.sample(list(genome_dictionary.keys()), 5000) # Lists the keys in the dictionary and takes 1000 random ones
    # the keys of the genome_dictionary is the bacteria IDs so it takes out 1000 bacteria IDs

# Distributions for the random bacteria (creates a new dictionary)   
genomes_distributions = {bacteria_id: genome_dictionary[bacteria_id] for bacteria_id in random_bacteria_ids} # Uses random_bacteria_ids to find 
    # the right keys in the dictionary and takes out the distributions for them 


# COMBINE GENE AND GENOME DIST ---------------------
# Make gene and genome distrbutions into dataframe
gene_df = pd.DataFrame(list(gene_distribution.items()), columns=["kmer", f"{gene_name}"])
gene_df.set_index("kmer", inplace=True) # inplace=True - modifies the gene_df instead of creating a new one
genomes_df = pd.DataFrame.from_dict(genomes_distributions, orient="index").T
genomes_df.fillna(0, inplace=True) # fill the NaN values with 0

# Join gene with genomes
distributions_df = gene_df.join(genomes_df, how="outer") # .join - joins based on the indexes, "outer" - takes all the kmers from both df
distributions_df.fillna(0, inplace=True)


# EUCLIDEAN DISTANCE -------------------------------
# Calculate eucldian distances between gene and genomes distributions
gene_vector = distributions_df[f"{gene_name}"].values[:, None]
eu_distances = np.linalg.norm(distributions_df.iloc[:, 1:].values - gene_vector, axis=0)

# Store bacteria id with the euclidiean distances 
euclidean_df = pd.DataFrame({
    'Bacteria_ID': distributions_df.columns[1:],  # The genome names (exclude the first column)
    'Euclidean_distance': eu_distances
})
    # euclidean_df contains all genomes we have chosen


# MATCH STATUS (from taxonomy) ----------------------
# Load Taxonomy results
path = "/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_1.csv" # taxonomy contains matches
taxonomy_df = pd.read_csv(path, sep=",", header=None)
filtered_df = taxonomy_df[taxonomy_df.iloc[:,0] == gene_name]   # Takes out the information for specific gene
filtered_df.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    # filtered_df contains all rows for the specified gene_name

# Classify the random bacteria IDs as match or no match if the gene exists or not
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
#print(match, no_match)

euclidean_df['Match_status'] = match_column # Add column with the match status


# ADD PHYLUM ----------------------------------------
# for all genomes, even the ones without a match
# Read full taxonomy file
taxonomy_file = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_taxonomy_df = pd.read_csv(taxonomy_file, sep="\t", header=None) 
full_taxonomy_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

euclidean_df = euclidean_df.merge(full_taxonomy_df[["Bacteria_ID", "Phylum"]], on="Bacteria_ID", how="left")


# ONLY PHYLA WITH MOST GENOMES ----------------------
top_phyla = euclidean_df["Phylum"].value_counts().head(6)
#print(top_phyla)
euclidean_filtered_df = euclidean_df[euclidean_df["Phylum"].isin(top_phyla.index)]
#print(euclidean_filtered_df.head())
    # euclidean_filtered_df - contains only the bacteria from the top phyla



# BINS FOR THE HISTOGRAM ----------------------------
# Make sure the bins are same size for all subplots 
nr_bins = 20 # Number of bins/bars 
min_value = euclidean_filtered_df["Euclidean_distance"].min()
max_value = euclidean_filtered_df["Euclidean_distance"].max() + 0.001 # so all values fall inside the max_value
bin_edges = np.linspace(min_value, max_value, nr_bins + 1)


# DOWNSAMPLE NO_MATCH -------------------------------
downsample_factor = 0.2 # keep 20% of the no_match bacteria, from each bin

downsampled_no_matches = [] # will become a list of dataframes

for phylum, phylum_df in euclidean_filtered_df.groupby("Phylum"): # phylum - name of phylum, phylum_df - df that has only rows from that phylum
    
    # create a df for match and a df for no_match
    matches_phylum_df = phylum_df[phylum_df["Match_status"] == "Match"]
    no_matches_phylum_df = phylum_df[phylum_df["Match_status"] == "No_match"]  

    # each bacteria is assigned to a bin
    bin_nr = np.digitize(no_matches_phylum_df["Euclidean_distance"], bin_edges) - 1  # subtracting 1 to get index starting from 0
    no_matches_phylum_df["bin_nr"] = bin_nr

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


# HISTOGRAM -----------------------------------------
# Create histogram with stacked bars
# plt.figure(figsize=(8, 5)) - not needed since we use sns.FacetGrid later

g = sns.FacetGrid(downsampled_df, col="Phylum",  col_order = top_phyla.index, sharey=False,
                  col_wrap = 3, height=4, aspect=1.2)  # euclidean_filtered_df before
g.map_dataframe(sns.histplot, x="Euclidean_distance", hue="Match_status", hue_order=["No_match", "Match"], multiple="stack", bins=bin_edges) # orange = match
#g.add_legend(title="Match Status")
g.set_axis_labels("Euclidean Distance", "Number of Bacteria")
#g.set_titles(col_template="{col_name}")
for ax, phylum in zip(g.axes.flat, top_phyla.index):
    ax.set_title(f"{phylum} (n={top_phyla[phylum]})")

g.set(xlim=(min_value - 0.001, max_value + 0.001))  
       
plt.savefig('/home/enyaa/gene_genome/histogram4.png') 
plt.close()



'''
Skräp?
palette = sns.color_palette(["#1F77B4", "#FF7F0E"])  # orange for "Match" and blue for "No Match"
# Apply the palette to the plot
sns.set_palette(palette)

sns.histplot(data=euclidean_filtered_df, x='Euclidean_distance', hue='Phylum', multiple='dodge', bins=nr_bins, palette='Set2', shrink=0.8)
    # hue='Match_status: Categorizes the data by Match_status, using different colours for each category and creates legend for this
    # multiple='stack': Stacks the histograms for different Match_status values instead of overlapping
sns.histplot(data=euclidean_filtered_df, x='Euclidean_distance',hue='Match_status', multiple='stack',  bins=nr_bins, palette='dark:#5A9_r') 

plt.xlabel('Eucledian distance')
plt.ylabel('Number of genomes')
plt.title('Title')
'''