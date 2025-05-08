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
            
            # Calculate how well the sequences align with this offset
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
    
    # Find best match for first motif
    offset1, score1, substring1 = find_best_motif_match(sequence, motif1)
    
    # Initial positions based on first motif alignment
    sequence_positions = {}
    for i, char in enumerate(sequence):
        sequence_positions[i] = offset1 + i
    
    # Find best match for second motif within reasonable spacing ranges
    best_score2 = 0
    best_offset2 = 0
    best_position2 = 0
    
    # Try different possible positions for the second motif
    for spacing in range(max_spacing + 1):
        expected_pos2 = len(motif1) + spacing
        
        # Try aligning different parts of the sequence with motif2
        for seq_start in range(max(0, offset1), min(len(sequence), offset1 + len(sequence))):
            remaining_seq = sequence[seq_start:]
            if not remaining_seq:
                continue
                
            offset2, score2, substring2 = find_best_motif_match(remaining_seq, motif2)
            adjusted_offset = seq_start + offset2
            
            # Calculate how close this is to our expected spacing
            spacing_penalty = abs(adjusted_offset - (offset1 + expected_pos2)) / max(1, max_spacing)
            adjusted_score = score2 * (1 - 0.1 * spacing_penalty)  # Penalize by spacing deviation
            
            if adjusted_score > best_score2:
                best_score2 = adjusted_score
                best_offset2 = offset2
                best_position2 = seq_start
    
    # Only use the second motif if we found a good match
    total_score = score1
    if best_score2 > len(motif2) * 0.3:  # At least 30% match
        total_score += best_score2
    
    # Place aligned data in the respective indices
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
    
    # Sort results by alignment score (descending)
    results.sort(key=lambda x: x[1], reverse=True)
    
    # Extract dataframes
    dfs = [result[0] for result in results]
    
    # Combine all partial DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

def find_dimer_motifs(sequences, min_length=4, max_length=6, max_spacing=5):
    """
    Try to identify two motifs in a set of sequences.
    This is a simplified approach - for real applications, you might want to use
    more sophisticated motif finding algorithms.
    """
    # Convert sequences to a position frequency matrix
    seq_length = len(sequences[0])
    matrix = np.zeros((4, seq_length))  # A, C, G, T
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    for seq in sequences:
        for i, base in enumerate(seq):
            if base in base_to_idx:
                matrix[base_to_idx[base], i] += 1
    
    # Normalize by column
    for j in range(seq_length):
        col_sum = sum(matrix[:, j])
        if col_sum > 0:
            matrix[:, j] /= col_sum
    
    # Find regions of high information content
    information_content = np.zeros(seq_length)
    for j in range(seq_length):
        for i in range(4):
            p = matrix[i, j]
            if p > 0:
                information_content[j] -= p * np.log2(p)
    
    # Find potential motif regions (high information content)
    motif_regions = []
    current_region = []
    
    threshold = np.mean(information_content) - 0.5 * np.std(information_content)
    
    for i, ic in enumerate(information_content):
        if ic < threshold:  # Lower information content means more conservation
            current_region.append(i)
        else:
            if len(current_region) >= min_length:
                motif_regions.append(current_region)
            current_region = []
    
    if current_region and len(current_region) >= min_length:
        motif_regions.append(current_region)
    
    # Keep only the most significant regions that are separated enough
    if len(motif_regions) < 2:
        # Not enough motifs found
        return None, None
    
    # Get consensus sequences for the two most significant motifs
    motif_regions.sort(key=len, reverse=True)
    motif_regions = motif_regions[:2]
    
    # Sort by position
    motif_regions.sort(key=lambda x: x[0])
    
    # Generate consensus sequences
    motifs = []
    for region in motif_regions:
        consensus = ""
        for pos in region:
            idx = np.argmax(matrix[:, pos])
            base = "ACGT"[idx]
            consensus += base
        motifs.append(consensus)
    
    if len(motifs) == 2:
        return motifs[0], motifs[1]
    else:
        return None, None