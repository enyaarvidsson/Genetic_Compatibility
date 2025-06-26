# 6. Count kmer for genes

from Bio import SeqIO
import subprocess
import os

k = 5 
threads = 8
fasta_file = "/storage/enyaa/nucleotide_fasta_protein_homolog_model.fasta" #"/home/enyaa/gene_genome/genes10.fasta"  
temp_directory = "/storage/enyaa/REVISED/KMER/gene_remove_temp" # used as kmc database & temporary directory
output_directory = "/storage/enyaa/REVISED/KMER/gene_kmer"

# Count k-mers for 10k genes separately
for record in SeqIO.parse(fasta_file, "fasta"):
    gene = record.id.split(" ")[0] 
    gene_name = gene.split("|")[5] # takes gene name for each gene
    #print(gene_name)
    
    if "/" in gene_name:
        gene_name = gene_name.replace("/", "?") # ? = /
        #print(gene_name)
    
    temp_fasta = os.path.join(temp_directory, f"{gene_name}.fasta") #temporary FASTA file for the current gene being processed
    output_prefix = os.path.join(temp_directory, f"{gene_name}_kmc") #kmc-specific file
    output_txt = os.path.join(output_directory, f"{gene_name}_kmer_counts.txt") #textfile
     #os.path.join: joins path and file
     
    #print(temp_fasta)
    
    # Write to temporary fasta file
    with open(temp_fasta, "w") as f:
        f.write(f">{record.id}\n{str(record.seq)}\n")
    
    # Run kmc on temprary fasta file
    subprocess.run([
        "kmc", "-fa", "-k" + str(k), "-t" + str(threads), "-ci1", "-cs10000", temp_fasta, output_prefix, temp_directory],
                   check=True)
    # "-fa": for FASTA file, "-b": non-canonical counts, "ci1": lowest num of counts to include (=1)
    # "-cs10000": maxmium number of counts (=10000), 
    
    #Dump k-mer counts to text file
    subprocess.run(["kmc_dump", output_prefix, output_txt], check=True)
    
    os.remove(temp_fasta) # Clean temporary fasta file
    
    
    
    







'''
blast_results_file = '/storage/enyaa/blast_results_10k.txt'
db_path = '/storage/enyaa/blast_db/genome_10k_db'

temp_dir = "/storage/enyaa/kmc_temp"
k = 5
threads = 1

with open(blast_results_file, 'r') as file:
    for line in file:
        fields = line.strip().split('\t')
        subject_id = fields[1] # Förmodligen 
        start = int(fields[8])
        end = int(fields[9])
        is_reverse = fields[1].startswith('-') #FEL
        
        cmd = f"blastdbcmd -db {db_path} -entry {subject_id} -range {start}-{end}"
        result = subprocess.run(cmd, shell=True, capture_out=True, text=True)
        
        matched_sequence = result.stdout.strip()
        
        if is_reverse:
            matched_sequence = str(Seq(matched_sequence).reverse_complement())
        
        # Count k-mers 
        subprocess.run(["kmc", "-k" +str(k), "-t" + str(threads), "ci1", ])
        

# Read file with taxonomy results 
taxonomy_file = "/home/enyaa/gene_genome/taxonomy_results_6.csv"
taxonomy_results_df = pd.read_csv(taxonomy_file, sep="\t")

# Parse FASTA file 
genomes_fasta_file = "/storage/enyaa/genomes_10k.fasta"
sequences = {record.id: record.seq for record in SeqIO.parse(genomes_fasta_file, "fasta")}



extracted_gene_seqs = []
for index, row in taxonomy_results_df.iterrows():
        genome_id = row['bacteria_id']
        start_position = row['Subject start']
        end_position = row['end_position']

        
        
        
        seq = sequences.get

'''
