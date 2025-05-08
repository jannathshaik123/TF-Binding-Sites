import os
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
from dimer_alignment_utils import generate_sorted_substrings, align_dimer_sequences_parallel, find_dimer_motifs
from using_meme import meme_analysis
from fasta_converter import fasta_format_converter
from seed_finder import get_top_sequences, seed_meme_analysis, read_html_pwm, calculate_consensus

def read_bzip_sequences(file_path, score_col=3):
    """
    Read bZIP dataset sequences and scores from the provided file format.
    
    Parameters:
    - file_path: Path to the bZIP 10mer file
    - score_col: Which score column to use (3 for z-score, 4 for rank-based)
    
    Returns:
    - tuple of (forward_sequences, forward_scores, reverse_sequences, reverse_scores)
    """
    try:
        # Read the file
        df = pd.read_csv(file_path, delimiter='\s+', header=None)
        
        # Check if file has the expected format
        if df.shape[1] < 4:
            print(f"Warning: File {file_path} does not have the expected 4 columns")
            return [], [], [], []
        
        # Extract forward and reverse sequences and their scores
        forward_sequences = df[0].tolist()
        reverse_sequences = df[1].tolist()
        forward_scores = df[score_col-1].tolist()  # Adjust for 0-based indexing
        reverse_scores = df[score_col-1].tolist()  # Using same scores for both
        
        return forward_sequences, forward_scores, reverse_sequences, reverse_scores
    
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return [], [], [], []

def is_heterodimer(filename):
    """Check if a file represents a heterodimer based on naming convention"""
    # If filename has format bZIP1_bZIP2_... it's a heterodimer
    parts = os.path.basename(filename).split('_')
    return len(parts) > 1 and not parts[0].endswith('.10mer')

def extract_tf_names(filename):
    """Extract TF names from the filename"""
    base = os.path.basename(filename)
    if is_heterodimer(filename):
        # For heterodimers, format is bZIP1_bZIP2_...
        parts = base.split('_')
        return f"{parts[0]}_{parts[1]}"
    else:
        # For homodimers, format is bZIP_...
        return base.split('_')[0]

def write_sequences_to_fasta(sequences, scores, output_file):
    """Write sequences and scores to a FASTA file"""
    with open(output_file, 'w') as file:
        for i, (sequence, score) in enumerate(zip(sequences, scores), 1):
            file.write(f">sequence_{i}|score={score:.4f}\n{sequence}\n")

def process_bzip_dimer(file_path, output_dir=None, max_spacing=5, top_n=50):
    """
    Process bZIP dimer binding data
    
    Parameters:
    - file_path: Path to the bZIP 10mer file
    - output_dir: Directory to save output files (default: based on filename)
    - max_spacing: Maximum spacing between dimer motifs
    - top_n: Number of top sequences to use for motif discovery
    """
    try:
        # Extract TF name from filename
        tf_name = extract_tf_names(file_path)
        dimer_type = "heterodimer" if is_heterodimer(file_path) else "homodimer"
        
        # Create output directory
        if output_dir is None:
            output_dir = f"bZIP_analysis/{tf_name}"
        os.makedirs(output_dir, exist_ok=True)
        
        print(f"Processing {dimer_type} {tf_name} from {file_path}")
        
        # Read sequences and scores
        forward_seqs, forward_scores, reverse_seqs, reverse_scores = read_bzip_sequences(file_path)
        
        if not forward_seqs:
            print(f"No sequences found in {file_path}")
            return False
        
        # Process forward and reverse sequences
        for direction, sequences, scores in [
            ("forward", forward_seqs, forward_scores),
            ("reverse", reverse_seqs, reverse_scores)
        ]:
            direction_dir = os.path.join(output_dir, direction)
            os.makedirs(direction_dir, exist_ok=True)
            
            # Get top sequences based on scores
            combined = list(zip(sequences, scores))
            sorted_combined = sorted(combined, key=lambda x: x[1], reverse=True)
            top_sequences = sorted_combined[:top_n]
            
            # Write top sequences to FASTA
            top_fasta_file = os.path.join(direction_dir, f"{tf_name}_top_{top_n}.fasta")
            write_sequences_to_fasta([seq for seq, _ in top_sequences], 
                                     [score for _, score in top_sequences], 
                                     top_fasta_file)
            
            # Extract just the sequences for motif discovery
            top_sequence_list = [seq for seq, _ in top_sequences]
            
            # Find potential dimer motifs
            print(f"Identifying motifs for {tf_name} ({direction})")
            motif1, motif2 = find_dimer_motifs(top_sequence_list, max_spacing=max_spacing)
            
            if not motif1 or not motif2:
                # If automatic detection fails, use MEME
                print(f"Automatic motif detection failed, using MEME for {tf_name} ({direction})")
                seed_meme_analysis(top_fasta_file, 0, top_n)
                html_file = os.path.join(direction_dir, f"Meme_of_top_{top_n}_Seeds/meme.html")
                
                if os.path.exists(html_file):
                    pwm_section = read_html_pwm(html_file)
                    consensus = calculate_consensus(pwm_section)
                    
                    # Try to split consensus into two motifs
                    motif_length = len(consensus) // 2
                    motif1 = consensus[:motif_length]
                    motif2 = consensus[motif_length:]
                    print(f"Potential motifs from MEME: {motif1} and {motif2}")
                else:
                    print(f"MEME HTML file not found at {html_file}, using default motifs")
                    # Use default short motifs if all else fails
                    motif1 = top_sequence_list[0][:5] 
                    motif2 = top_sequence_list[0][-5:]
            
            print(f"Using motifs: {motif1} and {motif2} for {tf_name} ({direction})")
            
            # Align all sequences using our dimer alignment method
            print(f"Aligning sequences for {tf_name} ({direction})")
            df_aligned = align_dimer_sequences_parallel(sequences, motif1, motif2, max_spacing)
            
            # Save the alignment results
            output_csv_path = os.path.join(direction_dir, f"{tf_name}_dimer_alignment.csv")
            df_aligned.to_csv(output_csv_path, index=False)
            print(f"Alignment data written to {output_csv_path}")
            
            # Extract and format aligned sequences
            seq_length = len(sequences[0])
            lower_bound = -seq_length
            upper_bound = 2 * seq_length
            
            # Get all available columns
            available_cols = [col for col in df_aligned.columns if lower_bound <= col < upper_bound]
            extracted_data = df_aligned[available_cols]
            formatted_data = extracted_data.apply(lambda row: ''.join(row.astype(str)), axis=1)
            
            # Write formatted sequences to a text file
            output_file_path = os.path.join(direction_dir, f"{tf_name}_dimer_consensus.txt")
            with open(output_file_path, 'w') as file:
                file.write(f"{motif1}-{motif2}\n")  # Write motifs as reference
                for line in formatted_data:
                    file.write(line.replace('-', 'N') + '\n')
            
            print(f"Formatted data written to {output_file_path}")
            
            # Convert to FASTA format for MEME analysis
            output_fasta_path = os.path.join(direction_dir, f"{tf_name}_dimer_consensus.fasta")
            fasta_format_converter(output_file_path, output_fasta_path)
            
            # Run MEME analysis on the aligned sequences
            meme_analysis(output_fasta_path, seq_length)
        
        return True
    
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return False

def process_bzip_directory(directory_path, pattern=".10mer.txt"):
    """Process all bZIP files in a directory"""
    success_count = 0
    failure_count = 0
    
    for root, _, files in os.walk(directory_path):
        for file in files:
            if pattern in file:
                file_path = os.path.join(root, file)
                print(f"Found bZIP file: {file_path}")
                
                success = process_bzip_dimer(file_path)
                if success:
                    success_count += 1
                else:
                    failure_count += 1
    
    print(f"Processing complete. Successful: {success_count}, Failed: {failure_count}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Process bZIP dimer binding data")
    parser.add_argument("--file", help="Path to a single bZIP 10mer file")
    parser.add_argument("--dir", help="Directory containing bZIP 10mer files")
    parser.add_argument("--top", type=int, default=50, help="Number of top sequences to use")
    parser.add_argument("--spacing", type=int, default=5, help="Maximum spacing between motifs")
    
    args = parser.parse_args()
    
    if args.file:
        process_bzip_dimer(args.file, top_n=args.top, max_spacing=args.spacing)
    elif args.dir:
        process_bzip_directory(args.dir)
    else:
        print("Please provide either --file or --dir argument")
        parser.print_help()