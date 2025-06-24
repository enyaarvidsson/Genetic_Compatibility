# Count k-mers (5-mers) for the genes

from Bio import SeqIO
import subprocess
import os



k = 5 
threads = 8
fasta_file = "/storage/enyaa/nucleotide_fasta_protein_homolog_model.fasta"  
temp_directory = "/storage/enyaa/FINAL/KMER/gene_remove_temp" # used as kmc database & temporary directory
output_directory = "/storage/enyaa/FINAL/KMER/gene_kmer"


for record in SeqIO.parse(fasta_file, "fasta"):
    gene = record.id.split(" ")[0] 
    gene_name = gene.split("|")[5] # takes gene name for each gene
    
    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?") # ? = /
    
    temp_fasta = os.path.join(temp_directory, f"{gene_name}.fasta") #temporary FASTA file for the current gene being processed
    output_prefix = os.path.join(temp_directory, f"{gene_name}_kmc") #kmc-specific file
    output_txt = os.path.join(output_directory, f"{gene_name}_kmer_counts.txt") #textfile
     #os.path.join: joins path and file
    
    # Write to temporary fasta file
    with open(temp_fasta, "w") as f:
        f.write(f">{record.id}\n{str(record.seq)}\n")
    
    # Run kmc on temprary fasta file
    subprocess.run([
        "kmc", "-fa", "-k" + str(k), "-t" + str(threads), "-ci1", "-cs10000", temp_fasta, output_prefix, temp_directory],
                   check=True)
    # "-fa": for FASTA file, "-b": non-canonical counts, "ci1": lowest num of counts to include (=1)
    # "-cs10000": maxmium number of counts (=10000), 
    
    # Dump k-mer counts to text file
    subprocess.run(["kmc_dump", output_prefix, output_txt], check=True)
    
    os.remove(temp_fasta) # Clean temporary fasta file
    
print("Done counting kmers for genes!")
    
