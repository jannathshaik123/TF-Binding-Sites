import pandas as pd
import numpy as np
from tqdm import tqdm
from itertools import product

def generate_sorted_substrings(s, min_length=3, max_length=8):
    """
    Generate sorted substrings from a string with length constraints.
    
    Args:
        s: Input string
        min_length: Minimum substring length
        max_length: Maximum substring length
    
    Returns:
        List of substrings sorted by length (descending) and lexicographically
    """
    substrings = [s[i:j] for i in range(len(s)) 
                 for j in range(i + min_length, min(i + max_length + 1, len(s) + 1))]
    substrings_sorted = sorted(substrings, key=lambda x: (-len(x), x))
    unique_substrings = []
    for substring in substrings_sorted:
        if substring not in unique_substrings:
            unique_substrings.append(substring)
    return unique_substrings

def find_dimer_motifs(sequences, min_motif_length=3, max_motif_length=6, 
                      min_spacing=0, max_spacing=10, top_n=5):
    """
    Find potential dimer motifs from a set of sequences.
    
    Args:
        sequences: List of DNA sequences
        min_motif_length: Minimum length of each motif
        max_motif_length: Maximum length of each motif
        min_spacing: Minimum spacing between motifs
        max_spacing: Maximum spacing between motifs
        top_n: Number of top motif pairs to return
    
    Returns:
        List of tuples (motif1, motif2, spacing, score)
    """
    # Count occurrences of all potential motifs
    motif_counts = {}
    for seq in sequences:
        substrings = generate_sorted_substrings(seq, min_motif_length, max_motif_length)
        for ss in substrings:
            motif_counts[ss] = motif_counts.get(ss, 0) + 1
    
    # Get top motifs based on occurrence
    top_motifs = sorted(motif_counts.items(), key=lambda x: x[1], reverse=True)[:50]
    top_motifs = [m[0] for m in top_motifs]
    
    # Evaluate all possible motif pairs with variable spacing
    motif_pair_scores = []
    
    for motif1, motif2 in tqdm(product(top_motifs, top_motifs), 
                              desc="Evaluating motif pairs",
                              total=len(top_motifs)**2):
        if motif1 == motif2:
            continue
            
        for spacing in range(min_spacing, max_spacing + 1):
            score = evaluate_motif_pair(sequences, motif1, motif2, spacing)
            motif_pair_scores.append((motif1, motif2, spacing, score))
    
    # Return top scoring motif pairs
    top_pairs = sorted(motif_pair_scores, key=lambda x: x[3], reverse=True)[:top_n]
    return top_pairs

def evaluate_motif_pair(sequences, motif1, motif2, spacing):
    """
    Evaluate how well a motif pair explains the sequences.
    
    Args:
        sequences: List of sequences
        motif1: First motif
        motif2: Second motif
        spacing: Spacing between motifs
    
    Returns:
        Score representing how well this motif pair explains the sequences
    """
    pattern = motif1 + '.' * spacing + motif2
    matches = 0
    partial_matches = 0
    
    for seq in sequences:
        # Check for exact match with specified spacing
        if motif1 in seq and motif2 in seq:
            idx1 = seq.find(motif1)
            if idx1 != -1:
                expected_idx2 = idx1 + len(motif1) + spacing
                if expected_idx2 + len(motif2) <= len(seq):
                    if seq[expected_idx2:expected_idx2+len(motif2)] == motif2:
                        matches += 1
                        continue
            
            # Check for presence of both motifs regardless of spacing
            partial_matches += 0.5
    
    # Score is weighted combination of exact and partial matches
    score = matches + 0.2 * partial_matches
    return score

def align_dimer_sequence(sequence_info):
    """
    Align a sequence based on two motifs with spacing.
    
    Args:
        sequence_info: Tuple of (sequence, motif1, motif2, spacing)
    
    Returns:
        DataFrame with aligned sequence
    """
    sequence, motif1, motif2, spacing = sequence_info
    max_length = len(sequence)
    lower_bound = 1 - max_length
    upper_bound = 2 * max_length - 1
    columns = [i for i in range(lower_bound, upper_bound)]
    df = pd.DataFrame('-', index=[0], columns=columns)
    
    # Try to find motif1 and motif2 with the specified spacing
    idx1 = sequence.find(motif1)
    
    if idx1 != -1:
        expected_idx2 = idx1 + len(motif1) + spacing
        
        if expected_idx2 + len(motif2) <= len(sequence) and sequence[expected_idx2:expected_idx2+len(motif2)] == motif2:
            # Place sequence in the dataframe, centered on the motifs
            indices = range(0, len(sequence))
            for i, val in zip(indices, sequence):
                if i in df.columns:
                    df.at[0, i] = val
            return df, True
    
    # If exact match not found, try to find best approximate match
    best_offset = 0
    best_match_score = -1
    
    for offset in range(-max_length, max_length):
        curr_score = 0
        
        # Check if motif1 can be found with this offset
        if motif1 in sequence:
            idx1 = sequence.find(motif1)
            curr_score += 1
            
            # Check if motif2 can be found at the expected position
            expected_idx2 = idx1 + len(motif1) + spacing
            if expected_idx2 + len(motif2) <= len(sequence) and expected_idx2 >= 0:
                overlap = 0
                for i in range(len(motif2)):
                    if expected_idx2 + i < len(sequence):
                        if sequence[expected_idx2 + i] == motif2[i]:
                            overlap += 1
                curr_score += overlap / len(motif2)
        
        if curr_score > best_match_score:
            best_match_score = curr_score
            best_offset = offset
    
    # Place sequence in the dataframe with best offset
    indices = range(best_offset, best_offset + len(sequence))
    for i, val in zip(indices, sequence):
        if i in df.columns:
            df.at[0, i] = val
    
    return df, False

def align_sequences_to_dimer(sequences, motif1, motif2, spacing):
    """
    Align multiple sequences based on dimer motifs.
    
    Args:
        sequences: List of sequences
        motif1: First motif
        motif2: Second motif
        spacing: Spacing between motifs
    
    Returns:
        DataFrame with aligned sequences
    """
    tasks = [(seq, motif1, motif2, spacing) for seq in sequences]
    
    aligned_dfs = []
    exact_matches = 0
    
    for task in tqdm(tasks, desc="Aligning sequences"):
        df, is_exact = align_dimer_sequence(task)
        aligned_dfs.append(df)
        if is_exact:
            exact_matches += 1
    
    print(f"Exact matches: {exact_matches}/{len(sequences)} ({exact_matches/len(sequences)*100:.2f}%)")
    
    # Combine all partial DataFrames
    combined_df = pd.concat(aligned_dfs, ignore_index=True)
    return combined_df

def identify_binding_sites(aligned_df, min_frequency=0.5):
    """
    Identify columns in the alignment with high conservation.
    
    Args:
        aligned_df: DataFrame with aligned sequences
        min_frequency: Minimum frequency for considering a column conserved
    
    Returns:
        List of column positions with high conservation
    """
    conserved_columns = []
    
    for col in aligned_df.columns:
        column_values = aligned_df[col].dropna()
        if len(column_values) == 0:
            continue
            
        # Count nucleotides in this column
        counts = {}
        for val in column_values:
            if val != '-':
                counts[val] = counts.get(val, 0) + 1
                
        # Check if any nucleotide is present at high frequency
        for nt, count in counts.items():
            freq = count / len(column_values)
            if freq >= min_frequency:
                conserved_columns.append((col, nt, freq))
                break
    
    return sorted(conserved_columns, key=lambda x: x[0])