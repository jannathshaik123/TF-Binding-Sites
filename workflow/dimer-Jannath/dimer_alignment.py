import pandas as pd
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from Bio import motifs
from Bio.Seq import Seq
import os
from alignment_utils import generate_sorted_substrings
from seed_finder import read_html_pwm, calculate_consensus
import subprocess
from collections import Counter

def read_sequences_from_file(file_path, col_index=0):
    """Read sequences from a file."""
    with open(file_path, 'r') as file:
        next(file)  # Skip header
        sequences = []
        scores = []
        for line in file:
            columns = line.split()
            if col_index == 0:
                sequences.append(columns[col_index])
            else:
                sequences.append(columns[col_index][::-1])
            scores.append(float(columns[2]))
    return sequences, scores

def find_common_motifs(sequences, min_len=4, max_len=8, top_n=20):
    """Find the most common motifs in the given sequences."""
    all_motifs = []
    for seq in sequences:
        for i in range(len(seq) - min_len + 1):
            for j in range(min_len, min(max_len + 1, len(seq) - i + 1)):
                motif = seq[i:i+j]
                all_motifs.append(motif)
    
    motif_counts = Counter(all_motifs)
    return motif_counts.most_common(top_n)

def align_sequence_to_dimer(sequence, motif1, motif2, spacings):
    """
    Align a sequence to a dimer pattern with variable spacing between motifs.
    Returns the best alignment and score.
    """
    best_alignment = None
    best_score = -1
    best_spacing = -1
    best_offset = 0
    
    for spacing in spacings:
        # Create pattern with current spacing
        pattern_len = len(motif1) + spacing + len(motif2)
        
        if len(sequence) < pattern_len:
            continue
        
        for i in range(len(sequence) - pattern_len + 1):
            # Check first motif
            motif1_score = sum(1 for j in range(len(motif1)) 
                              if i+j < len(sequence) and sequence[i+j] == motif1[j])
            
            # Check second motif
            motif2_start = i + len(motif1) + spacing
            motif2_score = 0
            if motif2_start + len(motif2) <= len(sequence):
                motif2_score = sum(1 for j in range(len(motif2))
                                  if sequence[motif2_start+j] == motif2[j])
            
            total_score = motif1_score + motif2_score
            
            if total_score > best_score:
                best_score = total_score
                best_spacing = spacing
                best_offset = i
                best_alignment = sequence[i:i+pattern_len]
    
    return best_alignment, best_score, best_spacing, best_offset

def align_sequences_to_dimer_parallel(sequences, motif1, motif2, spacings):
    """
    Align all sequences to the dimer motifs in parallel and create a DataFrame.
    """
    # Prepare tasks for parallel processing
    tasks = [(seq, motif1, motif2, spacings) for seq in sequences]
    
    with Pool(processes=cpu_count()) as pool:
        results = list(tqdm(pool.starmap(align_sequence_to_dimer, tasks), 
                           total=len(tasks), desc="Aligning sequences"))
    
    alignments, scores, used_spacings, offsets = zip(*results)
    
    # Determine the maximum length for aligned sequences
    max_len = max(len(align) for align in alignments if align is not None)
    
    # Create DataFrame for visualization
    df = pd.DataFrame('-', index=range(len(sequences)), columns=range(-max_len, max_len))
    
    for i, (alignment, spacing, offset) in enumerate(zip(alignments, used_spacings, offsets)):
        if alignment is not None:
            motif1_end = len(motif1)
            motif2_start = len(motif1) + spacing
            
            # Place aligned data in the DataFrame
            for j, nucleotide in enumerate(alignment):
                col_index = j - offset
                df.at[i, col_index] = nucleotide
    
    return df, used_spacings

def write_aligned_sequences_to_fasta(sequences, motif1, motif2, spacings, output_file):
    """
    Write aligned dimer sequences to a FASTA file.
    """
    results = []
    for seq in sequences:
        alignment, score, spacing, offset = align_sequence_to_dimer(seq, motif1, motif2, spacings)
        if alignment:
            results.append((alignment, score, spacing))
    
    # Sort by score
    results.sort(key=lambda x: x[1], reverse=True)
    
    with open(output_file, 'w') as file:
        for i, (alignment, score, spacing) in enumerate(results):
            file.write(f">sequence_{i+1}_spacing_{spacing}_score_{score}\n{alignment}\n")

def run_meme_for_dimer(fasta_file, output_dir, min_motif_width=4, max_motif_width=8):
    """
    Run MEME to find two motifs in the aligned sequences.
    """
    input_abs_path = os.path.abspath(fasta_file)
    input_dir, input_basename = os.path.split(input_abs_path)
    
    print(f"Running MEME on {input_abs_path}")
    try:
        subprocess.run(["docker", "run", "--rm", 
            "-v", f"{input_dir}:/data",  
            "memesuite/memesuite:latest", 
            "meme", f"/data/{input_basename}", 
            "-dna", 
            "-o",
            "-nostatus",
            "-minw", str(min_motif_width),
            "-maxw", str(max_motif_width),
            "-nmotifs", "2",  # Look for two motifs
            "-mod", "zoops",  # Zero or one occurrence per sequence
            "-objfun", "classic", 
            "-revcomp", 
            "-markov_order", "0", 
            "-o", f"/data/{output_dir}"],
           check=True)
        
        return f"{input_dir}/{output_dir}/meme.html"
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e}")
        return None

def extract_dimer_motifs_from_meme(html_file):
    """
    Extract the two motifs identified by MEME from the HTML output.
    """
    pwm_sections = []
    
    for i in range(2):  # Extract two motifs
        pwm_section = read_html_pwm(html_file, motif_index=i)
        if pwm_section:
            consensus = calculate_consensus(pwm_section)
            pwm_sections.append((consensus, pwm_section))
    
    return pwm_sections

def identify_optimal_spacing(sequences, motif1, motif2, min_spacing=0, max_spacing=10):
    """
    Identify the optimal spacing between two motifs.
    """
    spacing_scores = {}
    
    for spacing in range(min_spacing, max_spacing + 1):
        total_score = 0
        for seq in sequences:
            _, score, _, _ = align_sequence_to_dimer(seq, motif1, motif2, [spacing])
            total_score += score
        
        spacing_scores[spacing] = total_score / len(sequences)
    
    # Find the spacing with the highest score
    optimal_spacing = max(spacing_scores.items(), key=lambda x: x[1])[0]
    return optimal_spacing, spacing_scores

def main(input_file, output_dir, col_index=0, min_spacing=0, max_spacing=10):
    """
    Main function to run the dimer motif analysis.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Read sequences
    sequences, scores = read_sequences_from_file(input_file, col_index)
    
    # Find common motifs in sequences
    common_motifs = find_common_motifs(sequences)
    
    # Try combinations of the top motifs as potential dimers
    top_motifs = [motif for motif, _ in common_motifs[:10]]
    
    best_dimer = None
    best_score = -1
    best_spacing = -1
    
    for i, motif1 in enumerate(top_motifs):
        for j, motif2 in enumerate(top_motifs):
            if i == j:  # Skip same motif
                continue
                
            # Find optimal spacing
            optimal_spacing, spacing_scores = identify_optimal_spacing(
                sequences, motif1, motif2, min_spacing, max_spacing
            )
            
            score = spacing_scores[optimal_spacing]
            if score > best_score:
                best_score = score
                best_dimer = (motif1, motif2)
                best_spacing = optimal_spacing
    
    if best_dimer:
        motif1, motif2 = best_dimer
        print(f"Best dimer motifs: {motif1} - {best_spacing}bp - {motif2}")
        
        # Align sequences to the dimer motifs
        spacings = list(range(best_spacing-2, best_spacing+3))  # Allow some flexibility
        aligned_df, used_spacings = align_sequences_to_dimer_parallel(
            sequences, motif1, motif2, spacings
        )
        
        # Write to CSV
        output_csv = os.path.join(output_dir, "dimer_alignment.csv")
        aligned_df.to_csv(output_csv)
        print(f"Alignment written to {output_csv}")
        
        # Write to FASTA
        output_fasta = os.path.join(output_dir, "dimer_aligned.fasta")
        write_aligned_sequences_to_fasta(sequences, motif1, motif2, spacings, output_fasta)
        print(f"FASTA file written to {output_fasta}")
        
        # Run MEME on aligned sequences
        meme_output_dir = "meme_dimer_output"
        html_file = run_meme_for_dimer(output_fasta, meme_output_dir)
        
        if html_file:
            # Extract refined motifs from MEME output
            dimer_motifs = extract_dimer_motifs_from_meme(html_file)
            
            if len(dimer_motifs) == 2:
                refined_motif1, pwm1 = dimer_motifs[0]
                refined_motif2, pwm2 = dimer_motifs[1]
                
                print(f"Refined motifs:")
                print(f"Motif 1: {refined_motif1}")
                print(f"Motif 2: {refined_motif2}")
                print(f"Spacing: {best_spacing}")
    else:
        print("No good dimer motifs found")

if __name__ == "__main__":
    # Example usage
    input_file = r"C:\Users\rabia\Documents\.alphadock\TF-Binding-Sites\TF_Dimer\FOS_CEBPE_20N_lib9.10mer.txt"
    output_dir = "dimer_analysis_results"
    min_spacing = 0
    max_spacing = 10
    
    main(input_file, output_dir, col_index=0, min_spacing=min_spacing, max_spacing=max_spacing)