# How many streptomyces are in actinomycetota?

import pandas as pd


# Load one match_info file to get the taxonomic info for the filtered bacteria in our dataset
info_df = pd.read_csv("/storage/jolunds/FINAL/MATCH_INFO/NDM-1_match_info.tsv", sep="\t")
    # 72690 rows (bacteria), columns: Bacteria_ID, Domain, Phylum, Class, Order, Family, Genus, Species, Match_info
 

count_a = info_df[info_df['Phylum'] == 'Actinomycetota'].shape[0]
count_s = info_df[info_df['Genus'] == 'Streptomyces'].shape[0]
print(count_a)
print(count_s)

print(count_s / count_a) # 21 %
