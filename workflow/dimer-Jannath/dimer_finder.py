import pandas as pd
import numpy as np
from itertools import product
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from collections import Counter
import os
import subprocess
from Bio import motifs
from Bio.Seq import Seq

def find_potential_dimers(sequences, min_motif_len=4, max_motif_len=6, min_spacing=0, max_spacing=10):
    """
    Find potential dimer motifs in the given sequences.
    
    Args:
        sequences: List of DNA sequences
        min_motif_len: Minimum length of each motif
        max_motif_len: Maximum length of each motif
        min_spacing: Minimum spacing between motifs
        max_spacing: Maximum spacing between motifs
        
    Returns:
        List of potential dimer motifs with their scores
    """
    # Extract all possible k-mers from sequences
    all_kmers = []
    for k in range(min_motif_len, max_motif_len + 1):
        for seq in sequences:
            for i in range(len(seq) - k + 1):
                all_kmers.append(seq[i:i+k])
    
    # Count frequency of each k-mer
    kmer_counts = Counter(all_kmers)
    
    # Get the most common k-mers as potential motif candidates
    top_kmers = [kmer for kmer, count in kmer_counts.most_common(50)]
    
    # Test all possible combinations of k-mers with variable spacing
    potential_dimers = []
    
    for motif1 in tqdm(top_kmers, desc="Testing motif combinations"):
        for motif2 in top_kmers:
            for spacing in range(min_spacing, max_spacing + 1):
                dimer_score = evaluate_dimer(motif1, motif2, spacing, sequences)
                potential_dimers.append((motif1, motif2, spacing, dimer_score))
    
    # Sort by score in descending order
    potential_dimers.sort(key=lambda x: x[3], reverse=True)
    
    return potential_dimers

def evaluate_dimer(motif1, motif2, spacing, sequences):
    """
    Evaluate how well a dimer motif fits the sequences.
    
    Args:
        motif1: First motif
        motif2: Second motif
        spacing: Spacing between motifs
        sequences: List of sequences to evaluate against
        
    Returns:
        Score representing how well the dimer fits the sequences
    """
    dimer_pattern = motif1 + '.' * spacing + motif2
    pattern_len = len(dimer_pattern)
    
    # Count occurrences in sequences
    total_matches = 0
    for seq in sequences:
        if len(seq) < pattern_len:
            continue
            
        for i in range(len(seq) - pattern_len + 1):
            window = seq[i:i+pattern_len]
            matches = True
            
            for j in range(pattern_len):
                if dimer_pattern[j] == '.':
                    continue  # Skip spacing positions
                elif dimer_pattern[j] != window[j]:
                    matches = False
                    break
            
            if matches:
                total_matches += 1
                break  # Count only one match per sequence
    
    # Return percentage of sequences that contain the dimer
    return total_matches / len(sequences)

def generate_pwm_from_motif(motif_sequences):
    """
    Generate a Position Weight Matrix from aligned motif sequences
    """
    if not motif_sequences:
        return None
    
    instances = [Seq(seq) for seq in motif_sequences]
    m = motifs.create(instances)
    return m.counts

def extract_motif_instances(sequences, motif1, motif2, spacing, max_mismatches=1):
    """
    Extract instances of the motifs from the sequences allowing for mismatches
    """
    motif1_instances = []
    motif2_instances = []
    
    dimer_pattern = motif1 + '.' * spacing + motif2
    pattern_len = len(dimer_pattern)
    
    for seq in sequences:
        if len(seq) < pattern_len:
            continue
            
        for i in range(len(seq) - pattern_len + 1):
            window = seq[i:i+pattern_len]
            
            # Check first motif with allowed mismatches
            mismatches1 = sum(c1 != c2 for c1, c2 in zip(motif1, window[:len(motif1)]))
            
            # Check second motif with allowed mismatches
            motif2_start = len(motif1) + spacing
            mismatches2 = sum(c1 != c2 for c1, c2 in zip(motif2, window[motif2_start:motif2_start+len(motif2)]))
            
            if mismatches1 <= max_mismatches and mismatches2 <= max_mismatches:
                motif1_instances.append(window[:len(motif1)])
                motif2_instances.append(window[motif2_start:motif2_start+len(motif2)])
                break  # Only count one instance per sequence
    
    return motif1_instances, motif2_instances

def refine_motifs(sequences, motif1, motif2, spacing):
    """
    Refine the motifs by extracting instances and rebuilding the consensus
    """
    # Extract instances with some allowed mismatches
    motif1_instances, motif2_instances = extract_motif_instances(sequences, motif1, motif2, spacing)
    
    # Generate PWMs
    pwm1 = generate_pwm_from_motif(motif1_instances)
    pwm2 = generate_pwm_from_motif(motif2_instances)
    
    # Generate refined consensus motifs
    refined_motif1 = calculate_consensus_from_counts(pwm1)
    refined_motif2 = calculate_consensus_from_counts(pwm2)
    
    return refined_motif1, refined_motif2, pwm1, pwm2

def calculate_consensus_from_counts(counts):
    """
    Calculate consensus sequence from nucleotide counts
    """
    if not counts:
        return ""
        
    consensus = ""
    for i in range(len(counts['A'])):
        counts_at_pos = {base: counts[base][i] for base in "ACGT"}
        max_base = max(counts_at_pos.items(), key=lambda x: x[1])[0]
        consensus += max_base
    
    return consensus

def align_sequence_to_dimer(sequence, motif1, motif2, spacing):
    """
    Align a sequence to the dimer motifs
    """
    dimer_pattern = motif1 + '.' * spacing + motif2
    pattern_len = len(dimer_pattern)
    
    if len(sequence) < pattern_len:
        return None, -1
    
    best_match_pos = -1
    best_match_score = -1
    
    for i in range(len(sequence) - pattern_len + 1):
        window = sequence[i:i+pattern_len]
        
        # Score for motif1
        motif1_score = sum(c1 == c2 for c1, c2 in zip(motif1, window[:len(motif1)]))
        
        # Score for motif2
        motif2_start = len(motif1) + spacing
        motif2_score = sum(c1 == c2 for c1, c2 in zip(motif2, window[motif2_start:motif2_start+len(motif2)]))
        
        # Combined score
        total_score = motif1_score + motif2_score
        
        if total_score > best_match_score:
            best_match_score = total_score
            best_match_pos = i
    
    if best_match_pos != -1:
        return sequence[best_match_pos:best_match_pos+pattern_len], best_match_pos
    
    return None, -1

def align_sequences_to_dimer(sequences, motif1, motif2, spacing):
    """
    Align all sequences to the dimer motifs and create a DataFrame representation
    """
    aligned_seqs = []
    positions = []
    
    for seq in sequences:
        aligned_seq, pos = align_sequence_to_dimer(seq, motif1, motif2, spacing)
        if aligned_seq:
            aligned_seqs.append(aligned_seq)
            positions.append(pos)
        else:
            aligned_seqs.append('-' * (len(motif1) + spacing + len(motif2)))
            positions.append(-1)
    
    # Create a dataframe for visualization
    columns = list(range(len(motif1) + spacing + len(motif2)))
    df = pd.DataFrame(index=range(len(sequences)), columns=columns)
    
    for i, (seq, pos) in enumerate(zip(aligned_seqs, positions)):
        if pos != -1:
            for j, nucleotide in enumerate(seq):
                df.at[i, j] = nucleotide
        else:
            for j in range(len(columns)):
                df.at[i, j] = '-'
    
    return df

def create_logo_from_aligned_sequences(aligned_sequences, output_file):
    """
    Create a sequence logo from aligned sequences using WebLogo
    """
    # Placeholder for WebLogo implementation
    # In practice, you would use the WebLogo library or MEME Suite
    pass

def write_dimer_fasta(sequences, motif1, motif2, spacing, output_file):
    """
    Write aligned dimer sequences to a FASTA file
    """
    with open(output_file, 'w') as file:
        for i, seq in enumerate(sequences):
            aligned_seq, pos = align_sequence_to_dimer(seq, motif1, motif2, spacing)
            if aligned_seq:
                file.write(f">sequence_{i+1}\n{aligned_seq}\n")

def analyze_dimer_with_meme(input_file, motif1_len, motif2_len, min_spacing, max_spacing):
    """
    Use MEME to analyze a set of sequences for two motifs with variable spacing
    """
    input_abs_path = os.path.abspath(input_file)
    input_dir, input_basename = os.path.split(input_abs_path)
    output_dirname = f"{input_basename.replace('.fasta', '')}_meme_output"
    
    try:
        subprocess.run(["docker", "run", "--rm", 
            "-v", f"{input_dir}:/data",  
            "memesuite/memesuite:latest", 
            "meme", f"/data/{input_basename}", 
            "-dna", 
            "-o",
            "-nostatus",
            "-nmotifs", "2",  # Look for two motifs
            "-minw", str(min(motif1_len, motif2_len)), 
            "-maxw", str(max(motif1_len, motif2_len)),
            "-minsites", "10",
            "-mod", "zoops",  # Zero or one occurrence per sequence 
            "-objfun", "classic", 
            "-revcomp", 
            "-markov_order", "0", 
            "-o", f"/data/{output_dirname}"],
           check=True)
        
        return f"{input_dir}/{output_dirname}/meme.html"
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e}")
        return None

if __name__ == "__main__":
    # Example usage
    sequences = [
        "ATGACGCGTCAT",
        "ATGACTCGTCAT",
        "ATGACCCGTCAT",
        "CTGACGGGGTCAT",
        "ATGACGCGTCATA",
        # Add more sequences as needed
    ]
    
    # Find potential dimer motifs
    dimers = find_potential_dimers(sequences, min_spacing=1, max_spacing=5)
    
    # Get the top dimer candidate
    top_motif1, top_motif2, top_spacing, top_score = dimers[0]
    print(f"Top dimer candidate: {top_motif1}-{top_spacing}bp-{top_motif2}, Score: {top_score}")
    
    # Refine the motifs
    refined_motif1, refined_motif2, pwm1, pwm2 = refine_motifs(sequences, top_motif1, top_motif2, top_spacing)
    print(f"Refined motifs: {refined_motif1}-{top_spacing}bp-{refined_motif2}")
    
    # Align sequences to the dimer
    aligned_df = align_sequences_to_dimer(sequences, refined_motif1, refined_motif2, top_spacing)
    print(aligned_df)
    
    # Write aligned sequences to FASTA
    write_dimer_fasta(sequences, refined_motif1, refined_motif2, top_spacing, "dimer_aligned.fasta")
    
    # Analyze with MEME (would require Docker setup as in your original code)
    # meme_output = analyze_dimer_with_meme("dimer_aligned.fasta", len(refined_motif1), len(refined_motif2), top_spacing-2, top_spacing+2)