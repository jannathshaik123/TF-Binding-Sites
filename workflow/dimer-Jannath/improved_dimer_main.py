import os
import pandas as pd
import subprocess
from multiprocessing import Pool, cpu_count
from alignment_utils import generate_sorted_substrings, align_sequence
from using_meme import meme_analysis
from fasta_converter import fasta_format_converter
from seed_finder import get_top_sequences, seed_meme_analysis, read_html_pwm, calculate_consensus
from seq_logo_gen import sequence_logo_generator
import re
import numpy as np
from tqdm import tqdm

def read_sequences_from_file(file_path, col_index):
    """
    Read sequences from a space-separated file
    
    Parameters:
    - file_path: Path to the file
    - col_index: Index of the column containing sequences
    
    Returns:
    - Lists of sequences and scores
    """
    with open(file_path, 'r') as file:
        next(file)  # Skip header line
        sequences = []
        scores = []
        for line in file:
            columns = line.split()
            if col_index == 0:
                sequences.append(columns[col_index])
            else:
                sequences.append(columns[col_index][::-1])  # Reverse if col_index != 0
            scores.append(float(columns[2]))
    return sequences, scores

def read_sequences_from_fasta_file(file_path):
    """
    Read sequences from a FASTA file
    
    Parameters:
    - file_path: Path to the FASTA file
    
    Returns:
    - List of sequences
    """
    sequences = []
    with open(file_path, 'r') as file:
        current_seq = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line
        if current_seq:  # Add the last sequence
            sequences.append(current_seq)
    return sequences

def write_sequences_to_fasta(sequences, output_file):
    """
    Write sequences to a FASTA file
    
    Parameters:
    - sequences: List of (sequence, score) tuples
    - output_file: Path to the output file
    """
    with open(output_file, 'w') as file:
        for i, (sequence, score) in enumerate(sequences, 1):
            file.write(f">sequence_{i}\n{sequence}\n")

def align_sequences_to_dimer_parallel(sequences, motif1, motif2, spacing):
    """
    Align multiple sequences to a dimer motif in parallel
    
    Parameters:
    - sequences: List of sequences to align
    - motif1: First motif
    - motif2: Second motif
    - spacing: Spacing between motifs
    
    Returns:
    - DataFrame with all sequences aligned
    """
    tasks = [(seq, motif1, motif2, spacing) for seq in sequences]
    
    with Pool(processes=cpu_count()) as pool:
        dfs = pool.map(align_sequence_to_dimer, tasks)

    # Combine all partial DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

def align_sequence_to_dimer(sequence_info):
    """
    Align a sequence to a dimer motif
    
    Parameters:
    - sequence_info: Tuple of (sequence, motif1, motif2, spacing)
    
    Returns:
    - DataFrame with the aligned sequence
    """
    sequence, motif1, motif2, spacing = sequence_info
    
    # Calculate the maximum possible length
    max_length = len(sequence) + len(motif1) + spacing + len(motif2)
    
    # Create a DataFrame for the alignment
    # The range is to allow for shifts in both directions
    lower_bound = -len(sequence)
    upper_bound = max_length
    columns = [i for i in range(lower_bound, upper_bound)]
    df = pd.DataFrame('-', index=[0], columns=columns)
    
    best_score = 0
    best_start_pos = None
    
    # Try all possible positions of the sequence relative to the motifs
    combined_motif = motif1 + 'N' * spacing + motif2
    for start_pos in range(-len(sequence) + 1, len(combined_motif)):
        score = 0
        matches = 0
        
        # Count matches between sequence and motif1
        for i in range(len(motif1)):
            seq_pos = start_pos + i
            if 0 <= seq_pos < len(sequence):
                if sequence[seq_pos] == motif1[i]:
                    score += 1
                    matches += 1
        
        # Count matches between sequence and motif2
        for i in range(len(motif2)):
            seq_pos = start_pos + len(motif1) + spacing + i
            if 0 <= seq_pos < len(sequence):
                if sequence[seq_pos] == motif2[i]:
                    score += 1
                    matches += 1
        
        # If this is the best score so far, update the best position
        if score > best_score:
            best_score = score
            best_start_pos = start_pos
    
    if best_start_pos is not None:
        # Place the sequence in the DataFrame at the best position
        for i, char in enumerate(sequence):
            pos = best_start_pos + i
            if pos in df.columns:
                df.at[0, pos] = char
    
    return df

def extract_motifs_from_meme(html_file):
    """
    Extract multiple motifs from MEME HTML output
    
    Parameters:
    - html_file: Path to the MEME HTML output file
    
    Returns:
    - List of motifs extracted
    """
    from seed_finder import read_html_pwm, calculate_consensus
    
    try:
        pwm_section = read_html_pwm(html_file)
        if not pwm_section:
            return []
            
        consensus = calculate_consensus(pwm_section)
        if consensus:
            return [consensus]
        return []
    except Exception as e:
        print(f"Error extracting motifs from MEME HTML: {e}")
        return []

def analyze_dimer_binding_with_meme(TF, input_file_path, col_index=0, top_n=20, output_dir=None):
    """
    Analyze dimer binding using MEME for motif discovery
    
    Parameters:
    - TF: Transcription factor name
    - input_file_path: Path to the input file
    - col_index: Index of the column containing sequences
    - top_n: Number of top sequences to analyze
    - output_dir: Directory to store output files (default is DimerBinding/TF/col_X)
    
    Returns:
    - Dictionary with analysis results
    """
    # Create output directory if not provided
    if not output_dir:
        output_dir = f'DimerBinding/{TF}/col_{col_index}'
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Analyzing dimer binding for {TF}...")
    
    # Read sequences and scores
    sequences, scores = read_sequences_from_file(input_file_path, col_index)
    
    # Get top sequences
    top_sequences = get_top_sequences(sequences, scores, top_n=top_n)
    
    # Write top sequences to FASTA file
    output_fasta_file = os.path.join(output_dir, f'{TF}_top_{top_n}_sequences.fasta')
    write_sequences_to_fasta(top_sequences, output_fasta_file)
    
    # Run MEME analysis on top sequences to find motifs
    print(f"Running MEME analysis to find motifs...")
    seed_meme_analysis(output_fasta_file, col_index, top_n)
    
    # Path to MEME HTML output
    html_file = os.path.join(output_dir, f'Meme_of_top_{top_n}_Seeds/meme.html')
    
    # First get motif from initial MEME run
    motifs = extract_motifs_from_meme(html_file)
    
    # Run a second MEME analysis specifically looking for 2 motifs
    print("Running second MEME analysis for multiple motifs...")
    try:
        # Custom MEME analysis for multiple motifs
        subprocess.run([
            "docker", "run", "--rm", 
            "-v", f"{os.path.abspath(output_dir)}:/data",  
            "memesuite/memesuite:latest", 
            "meme", f"/data/{os.path.basename(output_fasta_file)}", 
            "-dna", 
            "-nostatus",
            "-maxw", "6",  # Smaller width for finding parts of dimers
            "-minw", "4",  # Smaller width for finding parts of dimers
            "-nmotifs", "2",  # Look for 2 motifs
            "-mod", "zoops", 
            "-objfun", "classic", 
            "-revcomp", 
            "-markov_order", "0", 
            "-o", f"/data/MultiMotif_Meme_Analysis"
        ], check=True)
        
        # Path to new MEME HTML output
        multi_html_file = os.path.join(output_dir, "MultiMotif_Meme_Analysis/meme.html")
        second_motifs = extract_motifs_from_meme(multi_html_file)
        
        # Add any new motifs found
        for motif in second_motifs:
            if motif not in motifs:
                motifs.append(motif)
                
    except Exception as e:
        print(f"Error running multiple motif MEME analysis: {e}")
    
    # Ensure we have at least 2 motifs
    if len(motifs) < 2:
        print("Not enough motifs found. Using most conserved regions as motifs...")
        # Find the most conserved regions in the sequences as a fallback
        from collections import Counter
        
        # Extract all substrings of length 4-6 from the sequences
        all_substrings = []
        for seq, _ in top_sequences:
            for i in range(len(seq) - 3):
                for j in range(4, 7):
                    if i + j <= len(seq):
                        all_substrings.append(seq[i:i+j])
        
        # Count frequency of each substring
        substring_counts = Counter(all_substrings)
        
        # Get the two most common substrings that don't overlap
        most_common = substring_counts.most_common(10)
        for substring, _ in most_common:
            # Add the substring if it doesn't overlap with existing motifs
            overlap = False
            for existing_motif in motifs:
                if substring in existing_motif or existing_motif in substring:
                    overlap = True
                    break
            
            if not overlap:
                motifs.append(substring)
                if len(motifs) >= 2:
                    break
    
    # Make sure we have exactly 2 motifs
    motifs = motifs[:2]
    if len(motifs) < 2:
        # Duplicate the first motif if we only found one
        motifs = [motifs[0], motifs[0]] if motifs else ["GATA", "GATA"]  # Default fallback
    
    motif1, motif2 = motifs[0], motifs[1]
    print(f"Identified motifs: {motif1} and {motif2}")
    
    # Find the optimal spacing between motifs
    print(f"Finding optimal spacing between motifs: {motif1} and {motif2}...")
    
    best_spacing = 0
    best_score = 0
    
    # Try different spacings
    for spacing in range(11):  # Try spacing from 0 to 10
        total_score = 0
        
        # Score each sequence against the dimer motif with this spacing
        for seq, _ in top_sequences:
            # Search for motifs in the sequence with up to 1 mismatch allowed
            for i in range(len(seq) - len(motif1) + 1):
                m1_matches = sum(a == b for a, b in zip(seq[i:i+len(motif1)], motif1))
                if m1_matches >= len(motif1) - 1:  # Allow 1 mismatch
                    expected_pos2 = i + len(motif1) + spacing
                    if expected_pos2 + len(motif2) <= len(seq):
                        m2_matches = sum(a == b for a, b in zip(seq[expected_pos2:expected_pos2+len(motif2)], motif2))
                        if m2_matches >= len(motif2) - 1:  # Allow 1 mismatch
                            total_score += m1_matches + m2_matches
        
        if total_score > best_score:
            best_score = total_score
            best_spacing = spacing
    
    print(f"Best dimer motif: {motif1}-{best_spacing * 'N'}-{motif2}")
    
    # Align all sequences to the best dimer motif
    print("Aligning all sequences to the best dimer motif...")
    all_alignments = align_sequences_to_dimer_parallel(sequences, motif1, motif2, best_spacing)
    
    # Save alignments to CSV
    output_csv_path = os.path.join(output_dir, f'{TF}_dimer_alignments.csv')
    all_alignments.to_csv(output_csv_path, index=False)
    print(f"Alignments saved to {output_csv_path}")
    
    # Extract aligned sequences for output
    aligned_sequences = []
    for _, row in all_alignments.iterrows():
        # Filter out '-' characters and join the sequence
        aligned_seq = ''.join([char for char in row if char != '-'])
        aligned_sequences.append(aligned_seq)
    
    # Write aligned sequences to text file
    output_txt_path = os.path.join(output_dir, f'{TF}_aligned_sequences.txt')
    with open(output_txt_path, 'w') as file:
        file.write(f"Motif 1: {motif1}\n")
        file.write(f"Spacing: {best_spacing}\n")
        file.write(f"Motif 2: {motif2}\n\n")
        for seq in aligned_sequences:
            file.write(f"{seq}\n")
    print(f"Aligned sequences saved to {output_txt_path}")
    
    # Convert to FASTA format
    output_fasta_path = os.path.join(output_dir, f'{TF}_aligned_sequences.fasta')
    fasta_format_converter(output_txt_path, output_fasta_path)
    print(f"FASTA file created at {output_fasta_path}")
    
    # Generate sequence logo
    print("Generating sequence logo...")
    try:
        sequence_logo_generator(TF, output_txt_path, col_index)
        print("Sequence logo generated successfully")
    except Exception as e:
        print(f"Error generating sequence logo: {e}")
    
    # Run full MEME analysis on aligned sequences
    print("Running MEME analysis on aligned sequences...")
    length_of_sequence = len(sequences[0])
    try:
        meme_analysis(output_fasta_path, length_of_sequence)
        print("MEME analysis completed successfully")
    except Exception as e:
        print(f"Error running MEME analysis: {e}")
    
    return {
        'motif1': motif1,
        'motif2': motif2,
        'spacing': best_spacing,
        'aligned_sequences': aligned_sequences,
        'output_dir': output_dir
    }

if __name__ == "__main__":
    # Example usage
    TF = 'FOS_CEBPE'
    input_file_path = r'C:\Users\rabia\Documents\.alphadock\TF-Binding-Sites\TF_Dimer\FOS_CEBPE_20N_lib9.10mer.txt'  # Update this path
    
    # Analyze dimer binding using MEME
    results = analyze_dimer_binding_with_meme(
        TF=TF, 
        input_file_path=input_file_path,
        col_index=0,
        top_n=20
    )
    
    # Print summary
    print("\nSummary of dimer motif found:")
    print(f"Motif 1: {results['motif1']}")
    print(f"Spacing: {results['spacing']}")
    print(f"Motif 2: {results['motif2']}")