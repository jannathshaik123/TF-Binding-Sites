import os
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import argparse

# Import utilities for monomer processing
from alignment_utils import generate_sorted_substrings as generate_sorted_substrings_monomer
from alignment_utils import align_sequence

# # Import utilities for dimer processing
# from dimer_alignment_utils import generate_sorted_substrings as generate_sorted_substrings_dimer
# from dimer_alignment_utils import align_dimer_sequences_parallel


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

def read_bzip_sequences(file_path, score_col=3):
    """Read sequences from bZIP format file."""
    try:
        df = pd.read_csv(file_path, delimiter='\s+', header=None)
        if df.shape[1] < 4:
            print(f"Warning: File {file_path} does not have the expected 4 columns")
            return [], []
        
        sequences = df[0].tolist()
        scores = df[score_col-1].tolist()  # Adjust for 0-based indexing
        return sequences, scores
    
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return [], []

def write_sequences_to_fasta(sequences, scores, output_file):
    """Write sequences and scores to a FASTA file."""
    with open(output_file, 'w') as file:
        for i, (sequence, score) in enumerate(zip(sequences, scores), 1):
            file.write(f">sequence_{i}|score={score:.4f}\n{sequence}\n")

def align_sequences_parallel(sequences, reference_sequence):
    """Align sequences in parallel (monomer mode)."""
    substrings = generate_sorted_substrings_monomer(reference_sequence)
    tasks = [(seq, reference_sequence, substrings) for seq in sequences]
    with Pool(processes=cpu_count()) as pool:
        dfs = pool.map(align_sequence, tasks)
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

def is_heterodimer(filename):
    """Check if a file represents a heterodimer based on naming convention."""
    parts = os.path.basename(filename).split('_')
    return len(parts) > 1 and not parts[0].endswith('.10mer')

def extract_tf_names(filename):
    """Extract TF names from the filename."""
    base = os.path.basename(filename)
    if is_heterodimer(filename):
        parts = base.split('_')
        return f"{parts[0]}_{parts[1]}"
    else:
        return base.split('_')[0]

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
        # Extract TF name from filename
        tf_name = extract_tf_names(file_path)
        dimer_type = "heterodimer" if is_heterodimer(file_path) else "homodimer"
        if output_dir is None:
            output_dir = f"DimerBinding/{tf_name}"
        os.makedirs(output_dir, exist_ok=True)
        print(f"Processing {dimer_type} {tf_name} from {file_path}")
        
        sequences, scores = read_bzip_sequences(file_path)
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

# def process_dimer(file_path, top_n=50, max_spacing=5, output_dir=None):
#     """
#     Process bZIP dimer binding data
    
#     Parameters:
#     - file_path: Path to the bZIP 10mer file
#     - top_n: Number of top sequences to use for motif discovery
#     - max_spacing: Maximum spacing between dimer motifs
#     - output_dir: Directory to save output files (default: based on filename)
#     """
#     try:
#         # Extract TF name from filename
#         tf_name = extract_tf_names(file_path)
#         dimer_type = "heterodimer" if is_heterodimer(file_path) else "homodimer"
#         if output_dir is None:
#             output_dir = f"DimerBinding/{tf_name}"
#         os.makedirs(output_dir, exist_ok=True)
#         print(f"Processing {dimer_type} {tf_name} from {file_path}")
#         sequences, scores = read_bzip_sequences(file_path)
#         if not sequences:
#             print(f"No sequences found in {file_path}")
#             return False
        
#         combined = list(zip(sequences, scores))
#         sorted_combined = sorted(combined, key=lambda x: x[1], reverse=True)
#         top_sequences = sorted_combined[:top_n]
#         top_fasta_file = os.path.join(output_dir, f"{tf_name}_top_{top_n}.fasta")
#         write_sequences_to_fasta([seq for seq, _ in top_sequences], 
#                                 [score for _, score in top_sequences], 
#                                 top_fasta_file)
        
#         top_sequence_list = [seq for seq, _ in top_sequences]
#         print(f"Identifying motifs for {tf_name}")
#         # motif1, motif2 = find_dimer_motifs(top_sequence_list, max_spacing=max_spacing)
#         # if not motif1 or not motif2:
#         # print(f"Automatic motif detection failed, using MEME for {tf_name}")
#         seed_meme_analysis(top_fasta_file, 0, top_n)
#         html_file = os.path.join(output_dir, f"Meme_of_top_{top_n}_Seeds/meme.html")
        
#         if os.path.exists(html_file):
#             pwm_section = read_html_pwm(html_file)
#             consensus = calculate_consensus(pwm_section)
#             motif_length = len(consensus) // 2
#             motif1 = consensus[:motif_length]
#             motif2 = consensus[motif_length:]
#             print(f"Potential motifs from MEME: {motif1} and {motif2}")
#         else:
#             print(f"MEME HTML file not found at {html_file}, using default motifs")
#             motif1 = top_sequence_list[0][:5] 
#             motif2 = top_sequence_list[0][-5:]
#         print(f"Using motifs: {motif1} and {motif2} for {tf_name}")
        
#         # Align all sequences using our dimer alignment method
#         print(f"Aligning sequences for {tf_name}")
#         df_aligned = align_dimer_sequences_parallel(sequences, motif1, motif2, max_spacing)
#         output_csv_path = os.path.join(output_dir, f"{tf_name}_dimer_alignment.csv")
#         df_aligned.to_csv(output_csv_path, index=False)
#         print(f"Alignment data written to {output_csv_path}")
        
#         seq_length = len(sequences[0])
#         lower_bound = -seq_length
#         upper_bound = 2 * seq_length
        
#         available_cols = [col for col in df_aligned.columns if lower_bound <= col < upper_bound]
#         extracted_data = df_aligned[available_cols]
#         formatted_data = extracted_data.apply(lambda row: ''.join(row.astype(str)), axis=1)
#         output_file_path = os.path.join(output_dir, f"{tf_name}_dimer_consensus.txt")
#         with open(output_file_path, 'w') as file:
#             file.write(f"{motif1}-{motif2}\n")  # Write motifs as reference
#             for line in formatted_data:
#                 file.write(line.replace('-', 'N') + '\n')
#         print(f"Formatted data written to {output_file_path}")
#         output_fasta_path = os.path.join(output_dir, f"{tf_name}_dimer_consensus.fasta")
#         fasta_format_converter(output_file_path, output_fasta_path)
#         meme_analysis(output_fasta_path, seq_length)
        
#         return True
    
#     except Exception as e:
#         print(f"Error processing dimer {file_path}: {e}")
#         return False


def process_directory(directory_path, file_type="auto", pattern=None):
    """Process all files in a directory based on specified criteria."""
    success_count = 0
    failure_count = 0
    monomer_pattern = "_8mers_top_enrichment"
    dimer_pattern = ".10mer"
    if pattern:
        monomer_pattern = pattern if file_type == "monomer" else monomer_pattern
        dimer_pattern = pattern if file_type == "dimer" else dimer_pattern
    
    
    for root, _, files in os.walk(directory_path):
        
        for file in files:
            file_path = os.path.join(root, file)
            if not (file.endswith(".txt") or file.endswith(".csv") or file.endswith(".tsv")):
                continue
            is_monomer = monomer_pattern in file
            is_dimer = dimer_pattern in file
            if file_type == "auto":
                if is_monomer:
                    print(f"Found monomer file: {file_path}")
                    success = process_monomer(file_path)
                elif is_dimer:
                    print(f"Found dimer file: {file_path}")
                    success = process_dimer(file_path)
                else:
                    print(f"Skipping file with unknown format: {file_path}")
                    continue
            elif file_type == "monomer" and (pattern is None or is_monomer):
                print(f"Processing as monomer: {file_path}")
                success = process_monomer(file_path)
            elif file_type == "dimer" and (pattern is None or is_dimer):
                print(f"Processing as dimer: {file_path}")
                success = process_dimer(file_path)
            else:
                continue
                
            if success:
                success_count += 1
            else:
                failure_count += 1
    
    print(f"Processing complete. Successful: {success_count}, Failed: {failure_count}")

def process_tf_list(tf_list, file_type="monomer", col_index=0):
    """Process a list of transcription factors."""
    for tf in tf_list:
        try:
            if file_type == "monomer":
                input_file_path = os.path.join(os.path.dirname(os.getcwd()),f'TF Binding Factors\{tf}\{tf}_8mers_top_enrichment.txt')
                process_monomer(input_file_path, col_index=col_index)
            # # elif file_type == "dimer":
            #     input_file_path = os.path.join(os.path.dirname(os.getcwd()),f'TF_Dimer\{tf}.10mer.txt')
            #     process_dimer(input_file_path)
        except Exception as e:
            print(f"Error processing {tf}: {e}")
            continue

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process transcription factor binding data for monomers and dimers")
    parser.add_argument("--file", help="Path to a single sequence file")
    parser.add_argument("--dir", help="Directory containing sequence files")
    parser.add_argument("--type", choices=["auto", "monomer", "dimer"], default="auto", 
                        help="Type of sequence data (auto-detect, monomer, or dimer)")
    parser.add_argument("--col", type=int, default=0, help="Column index for sequence data (monomer only)")
    parser.add_argument("--top", type=int, default=None, 
                        help="Number of top sequences to use (default: 20 for monomer, 50 for dimer)")
    parser.add_argument("--spacing", type=int, default=5, help="Maximum spacing between motifs (dimer only)")
    parser.add_argument("--pattern", help="File pattern to match when scanning directories")
    parser.add_argument("--tf_list", nargs='+', help="List of TF names to process (requires appropriate directory structure)")
    
    args = parser.parse_args()
    
    # Set default top_n based on type
    if args.top is None:
        args.top = 20 if args.type == "monomer" else 50
    
    if args.file:
        if args.type == "auto":
            # Auto-detect file type
            if ".10mer" in args.file:
                print(f"Auto-detected dimer file: {args.file}")
                process_dimer(args.file, top_n=args.top, max_spacing=args.spacing)
            else:
                print(f"Processing as monomer file: {args.file}")
                process_monomer(args.file, col_index=args.col, top_n=args.top)
        elif args.type == "monomer":
            process_monomer(args.file, col_index=args.col, top_n=args.top)
        # elif args.type == "dimer":
        #     process_dimer(args.file, top_n=args.top, max_spacing=args.spacing)
    elif args.dir:
        process_directory(args.dir, file_type=args.type, pattern=args.pattern)
    elif args.tf_list:
        process_tf_list(args.tf_list, file_type=args.type, col_index=args.col)
    else:
        print("Please provide either --file, --dir, or --tf_list argument")
        parser.print_help()