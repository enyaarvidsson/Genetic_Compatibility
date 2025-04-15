import pandas as pd
import re
from collections import Counter
from Bio.Seq import Seq

def tRNA_score_one_sided(gene_dist, genome_dist):
    #codons =[] #kanske behövs vid dubbelsidigt
    codon_scores = []
    for gene_codon, gene_value in gene_dist.items(): 
        if gene_codon in genome_dist:
            difference = max(0, gene_value-genome_dist[gene_codon])
            codon_scores.append(0.1*difference)
        else:
            first_two = gene_codon[:2]
            third = gene_codon[2]

            # Check wobble position 
            matched = False
            for genome_codon, genome_value in genome_dist.items():
                if genome_codon[:2] == first_two:
                    if third == "T" and genome_codon[2] == "C":
                        difference = max(0, gene_value - genome_value)
                        codon_scores.append(0.2*difference)
                        matched = True
                        break
                    elif third == "G" and genome_codon[2] == "A":
                        difference = max(0, gene_value - genome_value)
                        codon_scores.append(0.2*difference)
                        matched = True
                        break
                    elif genome_codon[2] == "T" and third != "G":
                        difference = max(0, gene_value - genome_value)
                        codon_scores.append(0.9*difference)
                        matched = True
                        break
            if not matched:
                codon_scores.append(gene_value)

    return sum(codon_scores)

def tRNA_score_two_sided(gene_dist, genome_dist):
    codon_scores = []
    for gene_codon, gene_value in gene_dist.items(): 
        if gene_codon in genome_dist:
            difference = abs(gene_value-genome_dist[gene_codon])
            codon_scores.append(0.1*difference)
        else:
            first_two = gene_codon[:2]
            third = gene_codon[2]

            # Check wobble position 
            matched = False
            for genome_codon, genome_value in genome_dist.items():
                if genome_codon[:2] == first_two:
                    if third == "T" and genome_codon[2] == "C":
                        difference = abs(gene_value - genome_value)
                        codon_scores.append(0.2*difference)
                        matched = True
                        break
                    elif third == "G" and genome_codon[2] == "A":
                        difference = abs(gene_value - genome_value)
                        codon_scores.append(0.2*difference)
                        matched = True
                        break
                    elif genome_codon[2] == "T" and third != "G":
                        difference = abs(gene_value - genome_value)
                        codon_scores.append(0.9*difference)
                        matched = True
                        break
            if not matched:
                codon_scores.append(gene_value)
                
    for genome_codon, genome_value in genome_dist.items():
        if genome_codon not in gene_dist:
            codon_scores.append(0.1*genome_value)

    return sum(codon_scores)

# Load file with gene_names (sorted)
gene_name_file = "/storage/enyaa/REVISED/gene_names.txt"
gene_names_df = pd.read_csv(gene_name_file, header=None, names=["Gene_name"])

gene_codons_file = '/storage/enyaa/REVISED/tRNA/codons_genes.csv'
gene_codons_df = pd.read_csv(gene_codons_file)

file = "/storage/enyaa/REVISED/KMER/euclidean_split_genes_filtered/euclidean_df_tet(Q).pkl"
filtered_df = pd.read_pickle(file)
bacteria_ids = filtered_df["Bacteria_ID"].tolist()

path = "/storage/shared/data_for_master_students/enya_and_johanna/genome_full_lineage.tsv"
full_lineage_df = pd.read_csv(path, sep="\t", header=None)
full_lineage_df.columns = ["Bacteria_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
phylum_mapping = full_lineage_df[["Bacteria_ID", "Phylum", "Species"]] # only the bacteria_id and the respective phylum & species



for gene_name in gene_names_df["Gene_name"]:
    gene_df = gene_codons_df[gene_codons_df['Gene_name'] == gene_name].reset_index(drop=True)
    stop_codons = ["TAA", "TAG", "TGA"]

    # Remove rows where 'Codon' is a stop codon
    gene_df = gene_df[~gene_df['Codon'].isin(stop_codons)].reset_index(drop=True)
    gene_df["Frequency"] = (
            gene_df["Count"] / gene_df["Count"].sum()
        )
    
    # Convert to a dictionary with Codon as the key and Count as the value
    gene_dict = gene_df.set_index("Codon")["Frequency"].to_dict()
    
    tRNA_score_list = []
    for bacteria_id in bacteria_ids:
        tRNA_file = f"/storage/jolunds/REVISED/tRNA/{bacteria_id}_trnascan.txt"
        tRNA_df = pd.read_csv(tRNA_file, sep="\t", comment="#", skiprows=3, header=None)
        tRNA_df.columns = ["Sequence_name", "tRNA_nr", "Begin", "End", "tRNA_type", "Anticodon", "Intron_begin", "Intron_end", "Score", "Comment"]
        
        genome_df = tRNA_df["Anticodon"].value_counts().reset_index()
        genome_df.columns = ["Anticodon", "Count"]
        
        # Reverse complement 
        genome_df["Codon"] = genome_df["Anticodon"].apply(
            lambda anticodon: str(Seq(anticodon).reverse_complement())
        )
        
        # Distribution
        genome_df["Frequency"] = (
            genome_df["Count"] / genome_df["Count"].sum()
        )
        
        # Change to dictionary - Codon as the key and Count as the value
        genome_dict = genome_df.set_index("Codon")["Frequency"].to_dict() 

        # tRNA score
        tRNA_value_one_sided = tRNA_score_one_sided(gene_dict, genome_dict)
        tRNA_value_two_sided = tRNA_score_two_sided(gene_dict, genome_dict)
        
        # Append results
        tRNA_score_list.append({"Bacteria_ID": bacteria_id, "tRNA_score_one_sided": tRNA_value_one_sided, 
                                "tRNA_score_two_sided": tRNA_value_two_sided})

    gene_tRNA_df = pd.DataFrame(tRNA_score_list)
    gene_tRNA_df = gene_tRNA_df.merge(phylum_mapping, on='Bacteria_ID', how='left')
    print(gene_tRNA_df)

    # Save
    if "/" in gene_name:
        gene_name = gene_name.rename("/", "?")
        
    save_path = f"/storage/enyaa/REVISED/tRNA/tRNA_score/tRNA_score_{gene_name}.csv"
    gene_tRNA_df.to_csv(save_path, index=False)


        