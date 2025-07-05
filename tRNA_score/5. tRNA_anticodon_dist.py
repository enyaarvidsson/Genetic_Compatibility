import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import re
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# Load bacteria ids
file = "./FINAL/filtered_bacteria.csv"
bacteria_df = pd.read_csv(file)
bacteria_ids = bacteria_df[["Bacteria_ID"]]

all_counts = []

for bacteria_id in tqdm(bacteria_ids["Bacteria_ID"], desc="Processing genomes"):
    tRNA_file = f"./FINAL/tRNA_SCAN/{bacteria_id}_trnascan.txt"
    if not os.path.exists(tRNA_file):
        continue
    
    tRNA_df = pd.read_csv(
        tRNA_file, sep="\t", comment="#", skiprows=3, header=None,
        names=["Sequence_name", "tRNA_nr", "Begin", "End", "tRNA_type",
               "Anticodon", "Intron_begin", "Intron_end", "Score", "Comment"]
    )
    anticodon_counts = tRNA_df["Anticodon"].value_counts()
    anticodon_counts["Bacteria_ID"] = bacteria_id
    all_counts.append(anticodon_counts)

# Create DataFrame with anticodons as columns
anticodon_counts_df = pd.DataFrame(all_counts).fillna(0)

# Merge with bacteria_ids_phylum
anticodon_counts_df = bacteria_df.merge(anticodon_counts_df, on="Bacteria_ID", how="right")

# Save dataframe
save_path = "./FINAL/anticodon_counts.csv"
anticodon_counts_df.to_csv(save_path)

# Create barplots for each phylum
file = "./FINAL/anticodon_counts.csv"
anticodon_counts_df = pd.read_csv(file, index_col=0)

id_columns = ["Bacteria_ID", "Phylum"]

# Identify anticodon columns (exclude ID columns)
possible_anticodons = [col for col in anticodon_counts_df.columns if col not in id_columns]

# Keep only anticodons that consist of exactly 3 characters, all A/C/T/G
valid_anticodons = [ac for ac in possible_anticodons if re.fullmatch(r"[ACGT]{3}", ac)]
anticodon_counts_df = anticodon_counts_df[id_columns + valid_anticodons]

top_phyla = ["Actinomycetota", "Bacillota", "Bacteroidota", "Campylobacterota", "Cyanobacteriota", "Pseudomonadota"]
for phylum in top_phyla:
    # Filter for phylum
    phylum_df = anticodon_counts_df[anticodon_counts_df["Phylum"] == phylum]

    # Melt into long format
    melted_df = phylum_df.melt(id_vars=["Bacteria_ID", "Phylum"], 
                                var_name="Anticodon", 
                                value_name="Count")
    
    melted_df["Anticodon"] = melted_df["Anticodon"].str.replace("T", "U")

    ######
    # Create custom legend items
    bar_patch = mpatches.Patch(color='C0', label='Average anticodon count')  # 'C0' is seaborn's default blue
    error_line = mlines.Line2D([], [], color='black', linewidth=1.5, label='Standard deviation')
    #####

    # Plot mean with standard deviation as error bars
    plt.figure(figsize=(8, 4))
    sns.barplot(
        data=melted_df,
        x="Anticodon",
        y="Count",
        errorbar="sd",
        capsize=0.1,    
        err_kws ={'linewidth': 1.5, 'color': 'black'} 
    )

    plt.xticks(rotation=90)
    plt.ylabel("", fontsize=14)
    plt.xlabel("Anticodon", fontsize=14)
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    
    ####
    plt.legend(handles=[bar_patch, error_line], fontsize=10)
    ###
    
    plt.tight_layout()
    plt.savefig(f"./FINAL/barplot_{phylum}.png")