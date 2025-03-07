# 4. Mobility classification

import pandas as pd

taxonomy_df = pd.read_csv("/storage/enyaa/REVISED/TAXONOMY/taxonomy_results_1.csv", header=None) # create a pandas dataframe
taxonomy_df.columns = ["Gene_name", "Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"] # add column names
taxonomy_levels = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# Count how many unqiue for each taxonomy level
taxonomy_counts = taxonomy_df.groupby("Gene_name")[taxonomy_levels].nunique() # group by Gene ID and count the unique values for different taxonomy levels, and count how many bacteria the gene is found in

# Cut-off mobile genes: (creates Boolean series)
mobile_genes = taxonomy_counts["Order"] > 1 # Has to be in more than one order

# Cut-off non-mobile genes: (creates Boolean series)
not_mobile_genes = taxonomy_counts["Genus"] == 1 # Has to be in only one genus, but can be in multiple species  

# Create the mobility column for the dataframe below
mobility_column = []
mobile = 0
not_mobile = 0
unclear = 0
for gene in taxonomy_counts.index:
    if gene in mobile_genes and mobile_genes[gene]: # Checks that gene is in row index, and also True in the Boolean series
        mobility_column.append("Mobile")
        mobile += 1
    elif gene in not_mobile_genes and not_mobile_genes[gene]:
        mobility_column.append("Not_mobile")
        not_mobile += 1
    else:
        mobility_column.append("Unclear")
        unclear += 1

# create a dataframe that says if the genes are mobile or not, and in how many bacteria it is found in
mobility_df = pd.DataFrame({
    "Gene_name": taxonomy_counts.index,
    "Mobility": mobility_column,
    "Number_of_bacteria": taxonomy_counts["Bacteria_ID"]
})
print(mobility_df.head())


# Save to csv-file
mobility_df.to_csv("/storage/jolunds/mobility_classification_1.csv", index=False)
print("File 'mobility_classification_1.csv' created successfully!")

# How many are mobile, not mobile
print("Mobile:", mobile)
print("Not mobile:", not_mobile)
print("Unclear:", unclear)

