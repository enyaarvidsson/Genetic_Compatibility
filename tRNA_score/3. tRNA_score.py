import pandas as pd
import re
from collections import Counter
from Bio.Seq import Seq
import os
import time
from tqdm import tqdm

def tRNA_score(gene_dist, genome_dist): # One-sided
    codon_scores = []
    max_score = 0
    max_codon = None

    for gene_codon, gene_value in gene_dist.items(): 
        score = 0  # default score
        if gene_codon in genome_dist:
            difference = max(0, gene_value - genome_dist[gene_codon])
            score = 1 * (difference ** 2)
        else:
            first_two = gene_codon[:2]
            third = gene_codon[2]
            matched = False
            for genome_codon, genome_value in genome_dist.items():
                if genome_codon[:2] == first_two:
                    if third == "T" and genome_codon[2] == "C":
                        difference = max(0, gene_value - genome_value)
                        score = 2 * (difference ** 2)
                        matched = True
                        break
                    elif third == "G" and genome_codon[2] == "A":
                        difference = max(0, gene_value - genome_value)
                        score = 2 * (difference ** 2)
                        matched = True
                        break
                    elif genome_codon[2] == "T" and third != "G":
                        difference = max(0, gene_value - genome_value)
                        score = 4 * (difference ** 2)
                        matched = True
                        break
            if not matched:
                score = 5 * (gene_value ** 2)

        codon_scores.append(score)

        if score > max_score:
            max_score = score
            max_codon = gene_codon

    return sum(codon_scores), max_score, max_codon

start_time = time.time()

# Load gene names
with open("./FINAL/gene_names.txt", "r") as f:
    all_genes = [line.strip() for line in f]

# Load file with counted codons for each gene 
gene_codons_file = './FINAL/codons_genes.csv'
gene_codons_df = pd.read_csv(gene_codons_file)

# Remove stop codons 
stop_codons = {"TAA", "TAG", "TGA"}
gene_codons_df = gene_codons_df[~gene_codons_df['Codon'].isin(stop_codons)]

# Load filtered bacteria 
file = "./FINAL/filtered_bacteria.csv"
bacteria_df = pd.read_csv(file)
bacteria_ids = bacteria_df["Bacteria_ID"].tolist()

# Compute frequencies for each gene
gene_freq_dict = {}
for gene_name, group in gene_codons_df.groupby("Gene_name"):
    freq_series = group["Count"] / group["Count"].sum()
    freq_dict = dict(zip(group["Codon"], freq_series))
    gene_freq_dict[gene_name] = freq_dict
    
# Compute frequencies for each genome
bacteria_codon_freq = {}
for bacteria_id in bacteria_ids:
    tRNA_file = f"./FINAL/tRNA_SCAN/{bacteria_id}_trnascan.txt"
    if not os.path.exists(tRNA_file):
        continue 
    
    tRNA_df = pd.read_csv(
        tRNA_file, sep="\t", comment="#", skiprows=3, header=None,
        names=["Sequence_name", "tRNA_nr", "Begin", "End", "tRNA_type",
               "Anticodon", "Intron_begin", "Intron_end", "Score", "Comment"]
        )
    
    # Value counts and reverse complement
    anticodon_counts = tRNA_df["Anticodon"].value_counts()
    valid_anticodon_counts = anticodon_counts[anticodon_counts.index.str.fullmatch(r'[ACGT]+')] # Remove anticodons containing other letters than ACTG
    
    codons = valid_anticodon_counts.index.to_series().apply(lambda ac: str(Seq(ac).reverse_complement())) # Mapping
    frequencies = valid_anticodon_counts / valid_anticodon_counts.sum()
    codon_freq = dict(zip(codons, frequencies))
    bacteria_codon_freq[bacteria_id] = codon_freq
    

for gene_name in tqdm(all_genes, desc="Processing genes"): 
    gene_dict = gene_freq_dict.get(gene_name, {})
    
    tRNA_score_list = []
    for bacteria_id in bacteria_ids:
        if bacteria_id not in bacteria_codon_freq:
            continue
        genome_dict = bacteria_codon_freq[bacteria_id]
    
        # tRNA score
        tRNA_value, tRNA_max_score, tRNA_max_codon = tRNA_score(gene_dict, genome_dict)

        # Append results with additional information
        tRNA_score_list.append({
            "Bacteria_ID": bacteria_id,
            "tRNA_score": tRNA_value,
            "Worst_case_tRNA": tRNA_max_score,
            "Worst_codon": tRNA_max_codon
        })

    gene_tRNA_df = pd.DataFrame(tRNA_score_list)
    gene_tRNA_df = gene_tRNA_df.merge(bacteria_df, on='Bacteria_ID', how='left')
   
    # Add match status and taxonomy
    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")
    
    match_file = f"./FINAL/MATCH_INFO/{gene_name}_match_info.tsv"
    gene_match_df = pd.read_csv(match_file, sep="\t")
    gene_match_df = gene_match_df[["Bacteria_ID", "Match_status"]]
    
    gene_tRNA_df = gene_tRNA_df.merge(gene_match_df, on=["Bacteria_ID"])
    
    # Save  
    directory = os.path.join('.', 'FINAL', 'tRNA_SCORE')
    os.makedirs(directory, exist_ok=True) 
    save_path = f"./FINAL/tRNA_SCORE/{gene_name}_tRNA_score.csv"
    gene_tRNA_df.to_csv(save_path, index=False)

end_time = time.time()
total_time = (end_time - start_time)/60
print(f"{total_time} minutes")


# Two-sided:
'''def tRNA_score_two_sided(gene_dist, genome_dist):
    codon_scores = []
    for gene_codon, gene_value in gene_dist.items(): 
        if gene_codon in genome_dist:
            difference = abs(gene_value-genome_dist[gene_codon])
            codon_scores.append(1 * (difference ** 2))
        else:
            first_two = gene_codon[:2]
            third = gene_codon[2]

            # Check wobble position 
            matched = False
            for genome_codon, genome_value in genome_dist.items():
                if genome_codon[:2] == first_two:
                    if third == "T" and genome_codon[2] == "C":
                        difference = abs(gene_value - genome_value)
                        codon_scores.append(2 * (difference ** 2))
                        matched = True
                        break
                    elif third == "G" and genome_codon[2] == "A":
                        difference = abs(gene_value - genome_value)
                        codon_scores.append(2 * (difference ** 2))
                        matched = True
                        break
                    elif genome_codon[2] == "T" and third != "G":
                        difference = abs(gene_value - genome_value)
                        codon_scores.append(4 * (difference ** 2))
                        matched = True
                        break
            if not matched:
                codon_scores.append(5 * (gene_value ** 2))
                
    for genome_codon, genome_value in genome_dist.items():
        if genome_codon not in gene_dist:
            codon_scores.append(1 * (genome_value ** 2))

    return sum(codon_scores)'''

