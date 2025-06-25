# Count k-mers (5-mers) for the length-adjusted genes

from Bio import SeqIO
from tqdm import tqdm
import random
import subprocess
import os


k = 5 
threads = 8
fasta_file = "/storage/enyaa/nucleotide_fasta_protein_homolog_model.fasta"  
temp_directory = "/storage/enyaa/FINAL/KMER/FOR_GENE_LENGTH/gene_remove_temp" # used as kmc database & temporary directory
output_directory = "/storage/enyaa/FINAL/KMER/FOR_GENE_LENGTH/gene_kmer_500bp"


# GENE LENGTH ----------------------------
for record in tqdm(SeqIO.parse(fasta_file, "fasta"), desc="Processing genes"):
    header = record.description
    positions = header.split("|")[3] 
    start = int(positions.split("-")[0])
    end = int(positions.split("-")[1])
    gene_length = abs(end-start) 
    gene_name = header.split("|")[-1].split(" ")[0]

    # skip genes that are shorter than 500 bp
    if gene_length < 500:
        continue

    if gene_length == 500:
        # take the whole sequence
        random_start = 0
    else:
        # get start position if index starts with 0
        random.seed(42)
        n = gene_length - 499
        random_start = random.randint(1, n) - 1 

    sequence = str(record.seq)  
    subseq = sequence[random_start:random_start+500] # the string will have 500 bp

    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?")

    # KMC ---------------------------------
    # create temp fasta-file and send to kmc
    temp_fasta = os.path.join(temp_directory, f"{gene_name}.fasta") # temporary FASTA file for the current gene being processed
    output_prefix = os.path.join(temp_directory, f"{gene_name}_kmc") # kmc-specific file
    output_txt = os.path.join(output_directory, f"{gene_name}_kmer_counts.txt") # textfile
        # os.path.join: joins path and file
    
    # Write to temporary fasta file
    with open(temp_fasta, "w") as f:
        f.write(f">{record.id}_subseq\n{subseq}\n")
    
    # Run kmc on temprary fasta file
    subprocess.run([
        "kmc", "-fa", "-k" + str(k), "-t" + str(threads), "-ci1", "-cs10000", temp_fasta, output_prefix, temp_directory],
                   check=True)
    # "-fa": for FASTA file, "-b": non-canonical counts, "ci1": lowest num of counts to include (=1)
    # "-cs10000": maxmium number of counts (=10000), 
    
    # Dump k-mer counts to text file
    subprocess.run(["kmc_dump", output_prefix, output_txt], check=True)
    
    os.remove(temp_fasta) # Clean temporary fasta file

print("Done counting kmers for the length-adjusted genes!")
