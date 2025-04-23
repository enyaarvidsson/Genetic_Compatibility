
# Count how many of each codon for the genes
# Save in a df

import pandas as pd
import re
from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO


# Count how many times each codon appears
def count_codons(sequence):
    sequence = re.sub(r'^ATGCRYWSKMBDHVN]', '', sequence.upper()) # Includes other letters so the codons will be correct
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    return Counter(codon for codon in codons) #if len(codon) == 3)


fasta_file = "/storage/enyaa/nucleotide_fasta_protein_homolog_model.fasta"
rows = []

for record in SeqIO.parse(fasta_file, "fasta"):
    header = record.id
    gene_name = header.split('|')[-1]

    codon_counts = count_codons(str(record.seq))
    
    for codon, count in codon_counts.items():
        rows.append({"Gene_name": gene_name, "Codon": codon, "Count": count})

gene_codons_df = pd.DataFrame(rows)

# Save to csv-file
output_file = '/storage/jolunds/REVISED/tRNA/codons_genes.csv'
gene_codons_df.to_csv(output_file, index=False)


# Remove codons with letters that are not ACGT
file = '/storage/jolunds/REVISED/tRNA/codons_genes.csv'
gene_codons_df = pd.read_csv(file)

incorrect_codons = ['GGN', 'NGG', 'GGS', 'GYC', 'KGG', 'AAR', 'ATY', 'CGY',
 'TCK', 'ANC', 'ANA', 'NAA', 'GNT', 'CNT', 'NTC', 'CGN', 'NTG', 'NAC', 'AAK', 'CMG',
 'TTY', 'CTY', 'ACK', 'YTG', 'CTN']

gene_codons_df = gene_codons_df[~gene_codons_df["Codon"].isin(incorrect_codons)]

save_path = "/storage/jolunds/REVISED/tRNA/codons_genes_filtered.csv"
gene_codons_df.to_csv(save_path)

