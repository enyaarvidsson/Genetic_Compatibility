import pandas as pd
import re
from collections import Counter
from Bio.Seq import Seq

def tRNA_score(codon, codon_value, codon_dist):
    if codon in codon_dist:
        difference = max(0, codon_value-codon_dist[codon])
        return 0.1*difference
    
    first_two = codon[:2]
    third = codon[2]

     # Check 0.5 rules first
    for dist_codon, dist_value in codon_dist.items():
        if dist_codon[:2] == first_two:
            if third == "T" and dist_codon[2] == "C":
                difference = max(0.05, codon_value - dist_value)
                return 0.2*difference
            elif third == "G" and dist_codon[2] == "A":
                difference = max(0.05, codon_value - dist_value)
                return 0.2*difference

    # Check 0.1 rule
    for dist_codon, dist_value in codon_dist.items():
        if dist_codon[:2] == first_two and dist_codon[2] == "T" and third != "G":
            difference = max(0.05, codon_value - dist_value)
            return 0.9*difference

    return codon_value


def count_codons(sequence):
    # Remove everything that isn't a base
    sequence = re.sub(r'[^ATGC]', '', sequence.upper()) 
    # Find codons
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    
    # Count codons
    codons_count = Counter(codons)
    
    codon_df = pd.DataFrame(codons_count.items(), columns=["Codon", "Count"])

    # Return the counts of anticodons
    return codon_df

# E.coli -------------------
Ecoli_file = "/home/jolunds/newtest/ecoli_trnascan_results.txt"
Ecoli_tRNA_df = pd.read_csv(Ecoli_file, sep="\t", comment="#", skiprows=3, header=None)
Ecoli_tRNA_df.columns = ["Sequence_name", "tRNA_nr", "Begin", "End", "tRNA_type", "Anticodon", "Intron_begin", "Intron_end", "Score", "Comment"]
Ecoli_anticodon_counts_df = Ecoli_tRNA_df["Anticodon"].value_counts().reset_index()
Ecoli_anticodon_counts_df.columns = ["Anticodon", "Ecoli_count"] 

# Reverse complement anticodons to get codons
Ecoli_anticodon_counts_df["Codon"] = Ecoli_anticodon_counts_df["Anticodon"].apply(
    lambda anticodon: str(Seq(anticodon).reverse_complement())
)
Ecoli_anticodon_counts_df["Frequency"] = (
    Ecoli_anticodon_counts_df["Ecoli_count"] / Ecoli_anticodon_counts_df["Ecoli_count"].sum()
)
# Change tp dictionary
Ecoli_codon_dist = dict(zip(Ecoli_anticodon_counts_df["Codon"], Ecoli_anticodon_counts_df["Frequency"]))

# B.subtilis ------------------
Bsubtilis_file = "/storage/enyaa/trna_scan/b_subtilis_trnascan_results.txt"
Bsubtilis_tRNA_df = pd.read_csv(Bsubtilis_file, sep="\t", comment="#", skiprows=3, header=None)
Bsubtilis_tRNA_df.columns = ["Sequence_name", "tRNA_nr", "Begin", "End", "tRNA_type", "Anticodon", "Intron_begin", "Intron_end", "Score", "Comment"]
Bsubtilis_anticodon_counts_df = Bsubtilis_tRNA_df["Anticodon"].value_counts().reset_index()
Bsubtilis_anticodon_counts_df.columns = ["Anticodon", "Bsubtilis_count"] 

# Reverse complement anticodons to get codons
Bsubtilis_anticodon_counts_df["Codon"] = Bsubtilis_anticodon_counts_df["Anticodon"].apply(
    lambda anticodon: str(Seq(anticodon).reverse_complement())
)
Bsubtilis_anticodon_counts_df["Frequency"] = (
    Bsubtilis_anticodon_counts_df["Bsubtilis_count"] / Bsubtilis_anticodon_counts_df["Bsubtilis_count"].sum()
)
# Change tp dictionary
Bsubtilis_codon_dist = dict(zip(Bsubtilis_anticodon_counts_df["Codon"], Bsubtilis_anticodon_counts_df["Frequency"]))

# B.fragilis --------------------
Bfragilis_file = "/storage/enyaa/trna_scan/B_fragilis_trnascan_results.txt"
Bfragilis_tRNA_df = pd.read_csv(Bfragilis_file, sep="\t", comment="#", skiprows=3, header=None)
Bfragilis_tRNA_df.columns = ["Sequence_name", "tRNA_nr", "Begin", "End", "tRNA_type", "Anticodon", "Intron_begin", "Intron_end", "Score", "Comment"]
Bfragilis_anticodon_counts_df = Bfragilis_tRNA_df["Anticodon"].value_counts().reset_index()
Bfragilis_anticodon_counts_df.columns = ["Anticodon", "Bfragilis_count"] 

# Reverse complement anticodons to get codons
Bfragilis_anticodon_counts_df["Codon"] = Bfragilis_anticodon_counts_df["Anticodon"].apply(
    lambda anticodon: str(Seq(anticodon).reverse_complement())
)
Bfragilis_anticodon_counts_df["Frequency"] = (
    Bfragilis_anticodon_counts_df["Bfragilis_count"] / Bfragilis_anticodon_counts_df["Bfragilis_count"].sum()
)
# Change tp dictionary
Bfragilis_codon_dist = dict(zip(Bfragilis_anticodon_counts_df["Codon"], Bfragilis_anticodon_counts_df["Frequency"]))

# NDM-6 -------------------------------
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

NDM6_codons = count_codons(NDM6_seqeunce)
NDM6_codons["Frequency"] = NDM6_codons["Count"] / NDM6_codons["Count"].sum()

NDM6_scores = []
for _, row in NDM6_codons.iterrows():
    codon = row["Codon"]
    freq = row["Frequency"]
    ecoli_codon_score = tRNA_score(codon, freq, Ecoli_codon_dist)
    bsubtilis_codon_score = tRNA_score(codon, freq, Bsubtilis_codon_dist)
    bfragilis_codon_score = tRNA_score(codon, freq, Bfragilis_codon_dist)
    NDM6_scores.append((
        codon,
        ecoli_codon_score,
        bsubtilis_codon_score,
        bfragilis_codon_score
    ))

NDM6_score_df = pd.DataFrame(NDM6_scores, columns=["Codon", "Ecoli_score", "Bsubtilis_score", "Bfragilis_score"])

NDM6_total_scores = NDM6_score_df[["Ecoli_score", "Bsubtilis_score", "Bfragilis_score"]].sum()
print(NDM6_total_scores)

# tet(Q) -------------------
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

tetQ_codons = count_codons(tetQ_seq)
tetQ_codons['Frequency'] = tetQ_codons['Count'] / tetQ_codons["Count"].sum()
tetQ_scores = []
for _, row in tetQ_codons.iterrows():
    codon = row['Codon']
    freq = row['Frequency']
    ecoli_codon_score = tRNA_score(codon, freq, Ecoli_codon_dist)
    bsubtilis_codon_score = tRNA_score(codon, freq, Bsubtilis_codon_dist)
    bfragilis_codon_score = tRNA_score(codon, freq, Bfragilis_codon_dist)
    
    tetQ_scores.append((
        codon,
        ecoli_codon_score,
        bsubtilis_codon_score,
        bfragilis_codon_score
    ))

tetQ_score_df = pd.DataFrame(tetQ_scores, columns=["Codon", "Ecoli_score", "Bsubtilis_score", "Bfragilis_score"])
tetQ_total_scores = tetQ_score_df[["Ecoli_score", "Bsubtilis_score", "Bfragilis_score"]].sum()
print(tetQ_total_scores)

#score_tetQ = tetQ_score_df['Score'].sum()/len(tetQ_score_df)
#print(score_tetQ)