import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from dimer_alignment_utils import (
    generate_sorted_substrings,
    find_dimer_motifs,
    align_sequences_to_dimer,
    identify_binding_sites
)
from weblogo import *
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_sequences_from_file(file_path, col_index=0):
    """
    Read sequences from a text file.
    
    Args:
        file_path: Path to the file containing sequences
        col_index: Column index containing sequences (0 for forward, 1 for reverse)
    
    Returns:
        Tuple of (sequences, scores)
    """
    try:
        with open(file_path, 'r') as file:
            next(file)  # Skip header
            sequences = []
            scores = []
            for line in file:
                columns = line.split()
                if len(columns) >= 3:  # Ensure we have enough columns
                    if col_index == 0:
                        sequences.append(columns[col_index])
                    else:
                        sequences.append(columns[col_index][::-1])  # Reverse for reverse strand
                    scores.append(float(columns[2]))
        return sequences, scores
    except Exception as e:
        print(f"Error reading sequences: {e}")
        return [], []

def read_sequences_from_fasta(file_path):
    """
    Read sequences from a FASTA file.
    
    Args:
        file_path: Path to the FASTA file
    
    Returns:
        List of sequences
    """
    sequences = []
    try:
        with open(file_path, 'r') as file:
            for record in SeqIO.parse(file, "fasta"):
                sequences.append(str(record.seq))
        return sequences
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return []

def get_top_sequences(sequences, scores, top_n=20):
    """
    Get top sequences based on scores.
    
    Args:
        sequences: List of sequences
        scores: List of scores
        top_n: Number of top sequences to return
    
    Returns:
        List of tuples (sequence, score)
    """
    combined = list(zip(sequences, scores))
    sorted_combined = sorted(combined, key=lambda x: x[1], reverse=True)
    top_sequences = sorted_combined[:top_n]
    return top_sequences

def write_sequences_to_fasta(sequences, output_file):
    """
    Write sequences to a FASTA file.
    
    Args:
        sequences: List of sequences or tuples (sequence, score)
        output_file: Path to the output FASTA file
    """
    with open(output_file, 'w') as file:
        for i, item in enumerate(sequences, 1):
            if isinstance(item, tuple):
                sequence, score = item
                file.write(f">sequence_{i}_score_{score}\n{sequence}\n")
            else:
                file.write(f">sequence_{i}\n{item}\n")

def run_meme_analysis(input_file, output_dir, min_width=6, max_width=8, 
                      nmotifs=2, mod="zoops", revcomp=True):
    """
    Run MEME analysis on a FASTA file.
    
    Args:
        input_file: Path to the input FASTA file
        output_dir: Directory to store MEME output
        min_width: Minimum motif width
        max_width: Maximum motif width
        nmotifs: Number of motifs to find
        mod: MEME model (zoops, oops, anr)
        revcomp: Consider reverse complement sequences
    """
    print(f"Running MEME analysis on: {input_file}")
    input_abs_path = os.path.abspath(input_file)
    input_dir, input_basename = os.path.split(input_abs_path)
    
    cmd = ["docker", "run", "--rm", 
           "-v", f"{input_dir}:/data",  
           "memesuite/memesuite:latest", 
           "meme", f"/data/{input_basename}", 
           "-dna", 
           "-oc", f"/data/{output_dir}",
           "-nostatus",
           "-minw", str(min_width), 
           "-maxw", str(max_width), 
           "-nmotifs", str(nmotifs), 
           "-mod", mod, 
           "-objfun", "classic"]
    
    if revcomp:
        cmd.append("-revcomp")
    
    try:
        subprocess.run(cmd, check=True)
        print(f"MEME analysis completed successfully. Output in: {os.path.join(input_dir, output_dir)}")
    except subprocess.CalledProcessError as e:
        print(f"MEME analysis failed: {e}")

def create_sequence_logo(aligned_df, output_file, title="Sequence Logo"):
    """
    Create a sequence logo from aligned sequences.
    
    Args:
        aligned_df: DataFrame with aligned sequences
        output_file: Path to save the output image
        title: Title for the sequence logo
    """
    # Convert DataFrame to proper format for WebLogo
    from weblogo import LogoData, LogoOptions, LogoFormat, png_print_formatter, read_seq_data
    from weblogo import std_color_schemes
    from io import StringIO
    
    # Create FASTA-format data from the DataFrame
    fasta_data = StringIO()
    for i, row in aligned_df.iterrows():
        seq_str = ''.join([val if val != '-' else 'N' for val in row.values])
        fasta_data.write(f">seq_{i}\n{seq_str}\n")
    
    # Reset position to start of StringIO
    fasta_data.seek(0)
    
    # Parse using weblogo's seq_data reader
    seqs = read_seq_data(fasta_data)
    
    # Generate logo
    logodata = LogoData.from_seqs(seqs)
    logooptions = LogoOptions()
    logooptions.title = title
    logooptions.color_scheme = std_color_schemes['classic']
    logooptions.yaxis_scale = 2.0
    logoformat = LogoFormat(logodata, logooptions)
    
    # Save logo to file
    png_bytes = png_print_formatter(logodata, logoformat)
    with open(output_file, 'wb') as fout:
        fout.write(png_bytes)
    
    print(f"Sequence logo created: {output_file}")

def main(input_file, output_dir, spacing_range=(0, 10), top_sequences=20):
    """
    Main function to run the dimer motif analysis pipeline.
    
    Args:
        input_file: Path to input file with sequences
        output_dir: Directory to store output files
        spacing_range: Tuple of (min_spacing, max_spacing)
        top_sequences: Number of top sequences to analyze
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read sequences from file
    if input_file.endswith('.fasta'):
        sequences = read_sequences_from_fasta(input_file)
        scores = [1.0] * len(sequences)  # Assign equal scores if not available
    else:
        sequences, scores = read_sequences_from_file(input_file)
    
    print(f"Read {len(sequences)} sequences from {input_file}")
    
    # Get top sequences
    top_seqs = get_top_sequences(sequences, scores, top_n=top_sequences)
    top_seq_only = [seq for seq, _ in top_seqs]
    
    # Write top sequences to FASTA file
    top_seqs_fasta = os.path.join(output_dir, "top_sequences.fasta")
    write_sequences_to_fasta(top_seqs, top_seqs_fasta)
    
    # Run MEME analysis to find potential motifs
    meme_output_dir = "meme_output"
    run_meme_analysis(top_seqs_fasta, meme_output_dir, 
                     min_width=3, max_width=8, nmotifs=2)
    
    print("\nUsing custom algorithm to find dimer motifs with variable spacing:")
    
    # Find potential dimer motifs using custom algorithm
    min_spacing, max_spacing = spacing_range
    dimer_motifs = find_dimer_motifs(top_seq_only, 
                                    min_motif_length=3, 
                                    max_motif_length=8,
                                    min_spacing=min_spacing, 
                                    max_spacing=max_spacing,
                                    top_n=5)
    
    print("\nTop potential dimer motifs:")
    for i, (motif1, motif2, spacing, score) in enumerate(dimer_motifs, 1):
        print(f"{i}. Motif1: {motif1}, Motif2: {motif2}, Spacing: {spacing}, Score: {score:.2f}")
    
    # Select the best motif pair for alignment
    best_motif1, best_motif2, best_spacing, _ = dimer_motifs[0]
    
    print(f"\nAligning sequences using best motif pair:")
    print(f"Motif1: {best_motif1}, Motif2: {best_motif2}, Spacing: {best_spacing}")
    
    # Align sequences using the best motif pair
    aligned_df = align_sequences_to_dimer(sequences, best_motif1, best_motif2, best_spacing)
    
    # Save aligned sequences to CSV
    aligned_csv = os.path.join(output_dir, f"aligned_dimer_spacing_{best_spacing}.csv")
    aligned_df.to_csv(aligned_csv, index=False)
    print(f"Aligned sequences saved to: {aligned_csv}")
    
    # Identify binding sites (conserved columns)
    binding_sites = identify_binding_sites(aligned_df, min_frequency=0.6)
    print("\nPotential binding sites (conserved columns):")
    for pos, nt, freq in binding_sites:
        print(f"Position {pos}: {nt} ({freq*100:.1f}%)")
    
    # Create sequence logo
    logo_file = os.path.join(output_dir, f"dimer_logo_spacing_{best_spacing}.png")
    create_sequence_logo(aligned_df, logo_file, 
                        title=f"Dimer Motif ({best_motif1}-{best_spacing}bp-{best_motif2})")
    
    print("\nDimer motif analysis completed successfully!")
    
    # Return the identified motifs and spacing for reference
    return best_motif1, best_motif2, best_spacing

if __name__ == "__main__":
    # Example usage
    input_file = r"C:\Users\rabia\Documents\.alphadock\TF-Binding-Sites\TF_Dimer\FOS_CEBPE_20N_lib9.10mer.txt"  # Replace with your input file
    output_dir = "dimer_analysis_results"
    
    # Specify the range of spacing to consider between motifs
    spacing_range = (0, 10)  # Try spacings from 0 to 10 bp
    
    main(input_file, output_dir, spacing_range)