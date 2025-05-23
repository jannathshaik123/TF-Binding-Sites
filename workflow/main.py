import os
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import argparse

from alignment_utils import generate_sorted_substrings, align_sequence
from using_meme import meme_analysis
from fasta_converter import fasta_format_converter
from seed_finder import get_top_sequences, seed_meme_analysis, read_html_pwm, calculate_consensus, find_dimer_motifs

def read_sequences_from_fasta_file(file_path, col_index=0):
    """Read sequences from a FASTA file."""
    sequences = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith(">"):
                sequences.append(line.strip())
    return sequences

def read_sequences_from_file(file_path, col_index=0):
    """Read sequences and scores from a tabular file."""
    with open(file_path, 'r') as file:
        next(file)
        sequences = []
        scores = []
        for line in file:
            columns = line.split()
            if col_index == 0:
                sequences.append(columns[col_index])
            else:
                sequences.append(columns[col_index][::-1])  # Reverse for col_index != 0
            scores.append(float(columns[2]))
    return sequences, scores

def write_sequences_to_fasta(sequences, scores, output_file):
    """Write sequences and scores to a FASTA file."""
    with open(output_file, 'w') as file:
        for i, (sequence, score) in enumerate(zip(sequences, scores), 1):
            file.write(f">sequence_{i}|score={score:.4f}\n{sequence}\n")

def align_sequences_parallel(sequences, reference_sequence):
    """Align sequences in parallel (monomer mode)."""
    substrings = generate_sorted_substrings(reference_sequence)
    tasks = [(seq, reference_sequence, substrings) for seq in sequences]
    with Pool(processes=cpu_count()) as pool:
        dfs = pool.map(align_sequence, tasks)
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

def process_monomer(file_path, col_index=0, top_n=20, output_dir=None):
    """
    Process transcription factor monomer binding data.
    
    Parameters:
    - file_path: Path to the sequence file
    - col_index: Column index for sequence data (0 or 1)
    - top_n: Number of top sequences to use
    - output_dir: Directory to save output files (default: based on filename)
    """
    try:
        tf_name = os.path.basename(file_path).split('_')[0]
        if output_dir is None:
            output_dir = f"MonomerBinding/{tf_name}/col_{col_index}"
        os.makedirs(output_dir, exist_ok=True)
        print(f"Processing monomer {tf_name} from {file_path}")
        
        sequences, scores = read_sequences_from_file(file_path, col_index)
        if not sequences:
            print(f"No sequences found in {file_path}")
            return False
        top_sequences = get_top_sequences(sequences, scores, top_n=top_n)
        output_fasta_file = os.path.join(output_dir, f'{tf_name}_top_{top_n}_sequences.fasta')
        write_sequences_to_fasta([seq for seq, _ in top_sequences], 
                                [score for _, score in top_sequences], 
                                output_fasta_file)
        
        seed_meme_analysis(output_fasta_file, col_index, top_n)
        html_file = os.path.join(output_dir, f"Meme_of_top_{top_n}_Seeds/meme.html")
        if not os.path.exists(html_file):
            print(f"MEME HTML file not found at {html_file}")
            return False
        pwm_section = read_html_pwm(html_file)
        reference_sequence = calculate_consensus(pwm_section)
        if not reference_sequence:
            print(f"No consensus sequence found for {tf_name}")
            return False
        print(f"Reference sequence for {tf_name}: {reference_sequence}")
        length_of_sequence = len(sequences[0])
        df_aligned = align_sequences_parallel(sequences, reference_sequence)
        
        output_csv_path = os.path.join(output_dir, f'{tf_name}_consensus.csv')
        df_aligned.to_csv(output_csv_path, index=False)
        print(f"Data written to {output_csv_path}")
        
        columns_to_extract = [col for col in range(length_of_sequence) if col in df_aligned.columns]
        extracted_data = df_aligned[columns_to_extract]
        formatted_data = extracted_data.apply(lambda row: ''.join(row.astype(str)), axis=1)
        
        output_file_path = os.path.join(output_dir, f'{tf_name}_consensus.txt')
        with open(output_file_path, 'w') as file:
            file.write(reference_sequence + '\n')
            for line in formatted_data:
                file.write(line + '\n')
        print(f"Data written to {output_file_path}")
        
        output_fasta_path = os.path.join(output_dir, f'{tf_name}_consensus.fasta')
        fasta_format_converter(output_file_path, output_fasta_path)
        meme_analysis(output_fasta_path, length_of_sequence)
        
        return True
    
    except Exception as e:
        print(f"Error processing monomer {file_path}: {e}")
        return False

def process_directory(directory_path, file_type="auto", pattern=None):
    """Process all files in a directory based on specified criteria."""
    success_count = 0
    failure_count = 0
    
    for root, _, files in os.walk(directory_path):
        
        for file in files:
            file_path = os.path.join(root, file)
            if not (file.endswith(".txt") or file.endswith(".csv") or file.endswith(".tsv")):
                continue
            print(f"Found monomer file: {file_path}")
            success = process_monomer(file_path)
                
            if success:
                success_count += 1
            else:
                failure_count += 1
    
    print(f"Processing complete. Successful: {success_count}, Failed: {failure_count}")

def process_tf_list(tf_list, file_type="monomer", col_index=0):
    """Process a list of transcription factors."""
    for tf in tf_list:
        try:
            input_file_path = os.path.join(os.path.dirname(os.getcwd()),f'TF Binding Factors\{tf}\{tf}_8mers.txt')
            process_monomer(input_file_path, col_index=col_index)
        except Exception as e:
            print(f"Error processing {tf}: {e}")
            continue

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process transcription factor binding data for monomers")
    parser.add_argument("--file", help="Path to a single sequence file")
    parser.add_argument("--dir", help="Directory containing sequence files")
    parser.add_argument("--col", type=int, default=0, help="Column index for sequence data (monomer only)")
    parser.add_argument("--top", type=int, default=None, 
                        help="Number of top sequences to use (default: 20)")
    parser.add_argument("--pattern", help="File pattern to match when scanning directories")
    parser.add_argument("--tf_list", nargs='+', help="List of TF names to process (requires appropriate directory structure)")
    
    args = parser.parse_args()
    
    # Set default top_n based on type
    if args.top is None:
        args.top = 20 if args.type == "monomer" else 50
    
    if args.file:
        process_monomer(args.file, col_index=args.col, top_n=args.top)
    elif args.dir:
        process_directory(args.dir, file_type=args.type, pattern=args.pattern)
    elif args.tf_list:
        process_tf_list(args.tf_list, file_type=args.type, col_index=args.col)
    else:
        print("Please provide either --file, --dir, or --tf_list argument")
        parser.print_help()
