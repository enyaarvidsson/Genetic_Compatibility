import pandas as pd
from collections import defaultdict
import os
import re
import time


def overlap_percentage(start1, end1, start2, end2):
    """
    Calculate overlap percentage between two coordinate ranges, 
    handling cases where sequences are on different strands.
    
    If sequences are on different strands, return 0 to indicate they should both be kept.
    """
    strand1 = "Forward"
    strand2 = "Forward"
    
    # Ensure start < end for both sequences
    if start1 > end1:
        start1, end1 = end1, start1 # Swap to ensure start < end
        strand1 = "Reverse"
    if start2 > end2:
        start2, end2 = end2, start2  # Swap to ensure start < end
        strand2 = "Reverse"
    
    if strand1 != strand2:
        return 0
    
    # Calculate overlap region
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap_length = max(0, overlap_end - overlap_start)  # Ensure no negative overlap

    # Calculate lengths of each sequence
    length1 = end1 - start1
    length2 = end2 - start2
    
    # Calculate overlap percentage for each sequence
    overlap_percentage1 = overlap_length / length1
    overlap_percentage2 = overlap_length / length2
    
    # Return the maximum overlap percentage
    return min(overlap_percentage1, overlap_percentage2)


def filter_blast_results(file_path):
    """Filter BLAST results based on GCA matches and overlap criteria."""
    with open(file_path, 'r') as f:
        lines = [line.strip().split() for line in f.readlines()]

    # Convert to DataFrame for easier processing
    df = pd.DataFrame(lines, columns=[
        "query", "subject", "identity", "alignment_length", "mismatches", "gap_opens",
        "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"
    ])

    # Convert necessary columns to numeric
    df[["s_start", "s_end", "bit_score"]] = df[["s_start", "s_end", "bit_score"]].apply(pd.to_numeric)

    # Sort by bit_score (higher is better)
    df = df.sort_values(by="bit_score", ascending=False)

    filtered_results = []
    seen_regions = {}
    removed_count = 0
    for _, row in df.iterrows():
        bacteria_id = row["subject"]
        s_start, s_end = row["s_start"], row["s_end"]

        if bacteria_id not in seen_regions:
            seen_regions[bacteria_id] = [(s_start, s_end)]
            filtered_results.append(row)
        else:
            keep = True
            for start, end in seen_regions[bacteria_id]:
                if overlap_percentage(start, end, s_start, s_end) > 0.2:
                    keep = False
                    break
            if keep:
                seen_regions[bacteria_id].append((s_start, s_end))
                filtered_results.append(row)
            else:
                removed_count += 1
    
    # Convert back to DataFrame and save the results
    filtered_df = pd.DataFrame(filtered_results)
    path = f"./FINAL/BLAST/filtered_blast_results.tsv"
    filtered_df.to_csv(path, sep="\t", index=False)

    
start_time = time.time()

path = "./FINAL/BLAST/blast_results.txt"
filter_blast_results(path)

end_time = time.time()
total_time = (end_time - start_time)/60

print(f"File created in {total_time} minutes")

