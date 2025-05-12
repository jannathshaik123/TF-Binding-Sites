import pandas as pd
import numpy as np
from tqdm import tqdm

def generate_sorted_substrings(s):
    """Generate all possible substrings of a string, sorted by length (descending)."""
    substrings = [s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1)]
    substrings_sorted = sorted(substrings, key=lambda x: (-len(x), x))
    unique_substrings = []
    for substring in substrings_sorted:
        if substring not in unique_substrings:
            unique_substrings.append(substring)
    return unique_substrings

def find_best_motif_match(sequence, motif, min_match_length=3):
    """Find the best match for a motif in a sequence."""
    substrings = generate_sorted_substrings(motif)
    
    best_match_score = 0
    best_offset = 0
    best_substring = ""
    
    for substring in substrings:
        if len(substring) < min_match_length:
            break  # Skip very short substrings
            
        if substring in sequence:
            idx = sequence.find(substring)
            substring_idx_in_motif = motif.find(substring)
            offset = substring_idx_in_motif - idx
            aligned_motif = motif[max(0, offset):]
            aligned_sequence = sequence[max(0, -offset):]
            
            min_length = min(len(aligned_motif), len(aligned_sequence))
            match_score = sum(1 for j in range(min_length) if aligned_motif[j] == aligned_sequence[j])
            
            if match_score > best_match_score:
                best_match_score = match_score
                best_offset = offset
                best_substring = substring
    
    return best_offset, best_match_score, best_substring

def align_dimer_sequence(sequence_info):
    """Align a sequence to two motifs (dimer) with variable spacing."""
    sequence, motif1, motif2, substrings1, substrings2, max_spacing = sequence_info
    
    max_length = len(sequence)
    lower_bound = -max_length
    upper_bound = 2 * max_length
    columns = [i for i in range(lower_bound, upper_bound)]
    df = pd.DataFrame('-', index=[0], columns=columns)
    

    offset1, score1, substring1 = find_best_motif_match(sequence, motif1)
    sequence_positions = {}
    for i, char in enumerate(sequence):
        sequence_positions[i] = offset1 + i
    best_score2 = 0
    best_offset2 = 0
    best_position2 = 0
    
    
    for spacing in range(max_spacing + 1):
        expected_pos2 = len(motif1) + spacing
        for seq_start in range(max(0, offset1), min(len(sequence), offset1 + len(sequence))):
            remaining_seq = sequence[seq_start:]
            if not remaining_seq:
                continue
                
            offset2, score2, substring2 = find_best_motif_match(remaining_seq, motif2)
            adjusted_offset = seq_start + offset2
            spacing_penalty = abs(adjusted_offset - (offset1 + expected_pos2)) / max(1, max_spacing)
            adjusted_score = score2 * (1 - 0.1 * spacing_penalty)  # Penalize by spacing deviation
            
            if adjusted_score > best_score2:
                best_score2 = adjusted_score
                best_offset2 = offset2
                best_position2 = seq_start
    
    total_score = score1
    if best_score2 > len(motif2) * 0.3:  # At least 30% match
        total_score += best_score2
    for i, val in enumerate(sequence):
        pos = sequence_positions.get(i, offset1 + i)
        if pos in df.columns:
            df.at[0, pos] = val
    
    return df, total_score

def align_dimer_sequences_parallel(sequences, motif1, motif2, max_spacing=5):
    """Align multiple sequences to a dimer motif pattern in parallel."""
    substrings1 = generate_sorted_substrings(motif1)
    substrings2 = generate_sorted_substrings(motif2)
    
    tasks = [(seq, motif1, motif2, substrings1, substrings2, max_spacing) for seq in sequences]
    
    results = []
    for task in tqdm(tasks, desc="Aligning sequences"):
        df, score = align_dimer_sequence(task)
        results.append((df, score))
    results.sort(key=lambda x: x[1], reverse=True)
    dfs = [result[0] for result in results]
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df
