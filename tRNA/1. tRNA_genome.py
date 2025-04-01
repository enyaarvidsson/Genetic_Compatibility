import subprocess
import pandas as pd

# Input FASTA file containing genome sequences
input_fasta = "/home/enyaa/gene_genome/B_fragilis.fna"
output_file = "/home/enyaa/gene_genome/B_fragilis_trnascan_results.txt"

'''# Run tRNAscan-SE
command = ["tRNAscan-SE", "-B", "-o", output_file, input_fasta] #-N output corresponding codon instead of anti-codon

try:
    subprocess.run(command, check=True)
    print(f"tRNAscan-SE completed successfully. Results saved to {output_file}")
except subprocess.CalledProcessError as e:
    print(f"Error running tRNAscan-SE: {e}")
 '''   

Ecoli_file = "/home/enyaa/gene_genome/e_coli_trnascan_results.txt"
Ecoli_tRNA_df = pd.read_csv(Ecoli_file, sep="\t", comment="#", skiprows=3, header=None)
Ecoli_tRNA_df.columns = ["Sequence_name", "tRNA_nr", "Begin", "End", "tRNA_type", "Anticodon", "Intron_begin", "Intron_end", "Score", "Comment"]

Ecoli_anticodon_counts_df = Ecoli_tRNA_df["Anticodon"].value_counts().reset_index()
Ecoli_anticodon_counts_df.columns = ["Anticodon", "Ecoli_count"] 


Bsubtilis_file = "/home/enyaa/gene_genome/b_subtilis_trnascan_results.txt"
Bsubtilis_tRNA_df = pd.read_csv(Bsubtilis_file, sep="\t", comment="#", skiprows=3, header=None)
Bsubtilis_tRNA_df.columns = ["Sequence_name", "tRNA_nr", "Begin", "End", "tRNA_type", "Anticodon", "Intron_begin", "Intron_end", "Score", "Comment"]

Bsubtilis_anticodon_counts_df = Bsubtilis_tRNA_df["Anticodon"].value_counts().reset_index()
Bsubtilis_anticodon_counts_df.columns = ["Anticodon", "Bsubtilis_count"] 

Bfragilis_file = "/home/enyaa/gene_genome/B_fragilis_trnascan_results.txt"
Bfragilis_tRNA_df = pd.read_csv(Bfragilis_file, sep="\t", comment="#", skiprows=3, header=None)
Bfragilis_tRNA_df.columns = ["Sequence_name", "tRNA_nr", "Begin", "End", "tRNA_type", "Anticodon", "Intron_begin", "Intron_end", "Score", "Comment"]

Bfragilis_anticodon_counts_df = Bfragilis_tRNA_df["Anticodon"].value_counts().reset_index()
Bfragilis_anticodon_counts_df.columns = ["Anticodon", "Bfragilis_count"]

#### GENE ---------
from collections import Counter
from Bio.Seq import Seq
import re
NDM6_seqeunce = """ATGGAATTGCCCAATATTATGCACCCGGTCGCGAAGCTGAGCACCGCATTAGCCGCTGCA
TTGATGCTGAGCGGGTGCATGCCCGGTGAAATCCGCCCGACGATTGGCCAGCAAATGGAA
ACTGGCGACCAACGGTTTGGCGATCTGGTTTTCCGCCAGCTCGCACCGAATGTCTGGCAG
CACACTTCCTATCTCGACATGCCGGGTTTCGGGGCAGTCGCTTCCAACGGTTTGATCGTC
AGGGATGGCGGCCGCGTGCTGGTGGTCGATACCGCCTGGACCGATGACCAGACCGCCCAG
ATCCTCAACTGGATCAAGCAGGAGATCAACCTGCCGGTCGCGCTGGCGGTGGTGACTCAC
GCGCATCAGGACAAGATGGGCGGTATGGACGCGCTGCATGCGGCGGGGATTGCGACTTAT
GCCAATGCGTTGTCGAACCAGCTTGCCCCGCAAGAGGGGATGGTTGCGGCGCAACACAGC
CTGACTTTCGCCGCCAATGGCTGGGTCGAACCAGCAACCGCGCCCAACTTTGGCCCGCTC
AAGGTATTTTACCCCGGCCCCGGCCACACCAGTGACAATATCACCGTTGGGATCGACGGC
ACCGACATCGCTTTTGGTGGCTGCCTGATCAAGGACAGCAAGGCCAAGTCGCTCGGCAAT
CTCGGTGATGCCGACACTGAGCACTACGCCGCGTCAGTGCGCGCGTTTGGTGCGGCGTTC
CCCAAGGCCAGCATGATCGTGATGAGCCATTCCGCCCCCGATAGCCGCGCCGCAATCACT
CATACGGCCCGCATGGCCGACAAGCTGCGCTGA"""

tetQ_seq = """GTGCGTTTCGACAATGCATCTATTGTAGTATATTATTGCTTAATCCAAATGAATATTATAAATTTAGGAATTCTTGCTCACATTGATGCA
GGAAAAACTTCCGTAACCGAGAATCTGCTGTTTGCCAGTGGAGCAACGGAAAAGTGCGGCCGTGTGGATAATGGTGACACCATAACAGAC
TCTATGGATATAGAGAAACGTAGAGGAATTACTGTTCGGGCTTCTACGACATCTATTATCTGGAATGGAGTGAAATGCAATATCATTGAC
ACTCCGGGACACATGGATTTTATTGCGGAAGTGGAGCGGACATTCAAAATGCTTGATGGAGCAGTCCTCATCTTATCCGCAAAGGAAGGC
ATACAAGCGCAAACAAAGTTGCTGTTCAATACTTTACAAAAACTGCAAATCCCGACAATTATATTTATCAATAAAATTGACCGTGACGGT
GTGAATTTAGAGCGTTTGTATCTGGATATAAAAACAAATCTGTCTCAAGATGTCCTGTTTATGCAAACTGTTGTCGATGGATTGGTTTAT
CCGATTTGCTCCCAAACATATATAAAGGAAGAATACAAAGAATTTGTATGCAACCATGACGACAATATATTAGAACGATATTTGGCGGAT
AGCGAAATTTCACCGGCTGATTATTGGAATACGATAATCGATCTTGTGGCAAAAGCCAAAGTCTATCCGGTACTACATGGATCAGCAATG
TTCAATATCGGTATCAATGAGTTGTTGGACGCCATCTCTTCTTTTATACTTCCTCCAGAATCAGTCTCAAACAGACTTTCAGCTTATCTC
TATAAGATAGAGCATGACCCCAAAGGACATAAAAGAAGTTTTCTAAAAATAATTGACGGAAGTCTGAGACTTCGAGACATTGTAAGAATC
AACGATTCGGAAAAATTCATCAAGATTAAAAATCTAAAGACTATTTATCAGGGCAGAGAGATAAATGTTGATGAAGTGGGGGCCAATGAT
ATCGCGATTGTAGAAGATATGGAAGATTTTCGAATCGGAGATTATTTAGGTACTAAACCTTGTTTGATTCAAGGGTTATCTCATCAGCAT
CCCGCTCTCAAATCCTCCGTCCGGCCAGACAGGTCCGAAGAGAGAAGCAAGGTGATATCCGCTCTGAATACATTGTGGATTGAAGACCCG
TCTTTGTCCTTTTCCATAAACTCATATAGTGATGAATTGGAAATCTCGTTATATGGTTTGACACAAAAGGAAATCATACAGACATTGCTG
GAAGAACGATTTTCCGTAAAGGTCCATTTTGATGAGATCAAGACTATCTACAAAGAACGACCTGTAAAAAAGGTCAATAAGATTATTCAG
ATCGAAGTGCCACCCAACCCTTACTGGGCCACAATAGGGCTGACGCTTGAACCCTTGCCGTTAGGGACAGGGTTGCAAATCGAAAGTGAC
ATCTCCTATGGTTATCTGAACCATTCTTTTCAAAATGCCGTTTTTGAAGGGATTCGTATGTCTTGCCAATCTGGTTTACATGGATGGGAA
GTGACTGATCTGAAAGTAACTTTTACTCAAGCCGAGTATTATAGCCCGGTAAGTACACCTGCTGATTTCAGACAGCTGACCCCTTATGTC
TTCAGGCTGGCCTTGCAACAGTCAGGTGTGGACATTCTCGAACCGATGCTCTATTTTGAGTTGCAGATACCCCAAGCGGCAAGTTCCAAA
GCTATTACAGATTTGCAAAAAATGATGTCTGAGATTGAAGACATCAGTTGCAATAATGAGTGGTGTCATATTAAAGGGAAAGTTCCATTA
AATACAAGTAAAGACTACGCCTCAGAAGTAAGTTCATACACTAAGGGCTTAGGCGTTTTTATGGTCAAGCCATGCGGGTATCAAATAACA
AAAGGCGATTATTCTGATAATATCCGCATGAACGAAAAAGATAAACTTTTATTCATGTTCCAAAAATCAATGTCATCAAAATAA"""

def count_anticodons(sequence):
    # Remove everything that isn't a base
    sequence = re.sub(r'[^ATGC]', '', sequence.upper()) 
    # Find codons
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    
    # Convert to anticodns using reverse complement
    anticodons = [str(Seq(codon).reverse_complement()) for codon in codons]
    
    # Count anticodons
    anticodon_counts = Counter(anticodons)
    
    anticodon_df = pd.DataFrame(anticodon_counts.items(), columns=["Anticodon", "Gene_count"])

    # Return the counts of anticodons
    return anticodon_df
    

NDM6_anticodons = count_anticodons(NDM6_seqeunce)
tetQ_anticodons = count_anticodons(tetQ_seq)
#print(NDM6_anticodons)

#compare_df = tetQ_anticodons.merge(Ecoli_anticodon_counts_df, on=['Anticodon'], how='outer').fillna(0)
compare_df = (
    NDM6_anticodons
    .merge(Ecoli_anticodon_counts_df, on="Anticodon", how="outer")
    .merge(Bsubtilis_anticodon_counts_df, on="Anticodon", how="outer")
    .merge(Bfragilis_anticodon_counts_df, on="Anticodon", how="outer")
    .fillna(0)
)
compare_df['Ecoli_count'] = compare_df['Ecoli_count'].astype(int)
compare_df['Bsubtilis_count'] = compare_df['Bsubtilis_count'].astype(int)
compare_df['Bfragilis_count'] = compare_df['Bfragilis_count'].astype(int)
compare_df['Gene_count'] = compare_df['Gene_count'].astype(int)

compare_df["Ecoli_compatibility"] = compare_df.apply(lambda row: max(0, row["Gene_count"] - row["Ecoli_count"]) if row["Ecoli_count"] == 0 else 0, axis=1)
compare_df["Bsubtilis_compatibility"] = compare_df.apply(lambda row: max(0, row["Gene_count"] - row["Bsubtilis_count"]) if row["Bsubtilis_count"] == 0 else 0, axis=1)
compare_df["Bfragilis_compatibility"] = compare_df.apply(lambda row: max(0, row["Gene_count"] - row["Bfragilis_count"]) if row["Bfragilis_count"] == 0 else 0, axis=1)

ecoli_score = compare_df["Ecoli_compatibility"].sum()/compare_df["Gene_count"].sum()
bsubtilis_score = compare_df["Bsubtilis_compatibility"].sum()/compare_df["Gene_count"].sum()
bfragilis_score = compare_df["Bfragilis_compatibility"].sum()/compare_df["Gene_count"].sum()

#print(compare_df.head(30))
#print(compare_df.tail(27))

print(ecoli_score)
print(bsubtilis_score)
print(bfragilis_score)
