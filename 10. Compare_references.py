import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

path_incomp = "/storage/enyaa/incompatible_reference.tsv"
path_comp = "/storage/enyaa/compatible_reference.tsv"

incomp_df = pd.read_csv(path_incomp, sep="\t")
comp_df = pd.read_csv(path_comp, sep="\t")

# Merge:
merged_df = pd.concat([incomp_df, comp_df], axis=0, ignore_index=True )
print(merged_df)

# Histogram:
# Create histogram with stacked bars
plt.figure(figsize=(8, 5))

sns.histplot(data=merged_df, x='Euclidean_distance', hue='Gene_name', multiple='stack',  bins=20) 

plt.xlabel('Euclidean distance')
plt.ylabel('Number of genomes')
plt.title('Title')

plt.savefig('/home/enyaa/gene_genome/histogram_merged.png') 
plt.close()