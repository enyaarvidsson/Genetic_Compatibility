# Take the gene names from the file "nucleotide_fasta_protein_homolog_model.fasta"
# and sort them alphabetically. 

import re
from Bio import SeqIO


# headers in the fasta file
filepath = "/storage/enyaa/nucleotide_fasta_protein_homolog_model.fasta"
fasta_headers = [record.description for record in SeqIO.parse(filepath, "fasta")]

# loop through each header to extract gene names
extracted_data = []
for header in fasta_headers:
    parts = header.split("|")
    gene_name = parts[-1].split("[")[0].strip() # get the gene name
    
    extracted_data.append(gene_name)

# function to sort the gene names
def natural_sort_key(s):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]

# sort the gene names using natural sorting
extracted_data.sort(key=natural_sort_key)

# save the sorted gene names in txt file
output_file = "/storage/enyaa/FINAL/gene_names.txt"
with open(output_file, "w") as f:
    f.write("\n".join(extracted_data))

print(f"Gene names have been sorted")

