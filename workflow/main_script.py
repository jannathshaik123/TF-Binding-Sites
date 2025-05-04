import pandas as pd
from multiprocessing import Pool, cpu_count
from  alignment_utils import generate_sorted_substrings, align_sequence
from  using_meme import meme_analysis
from  fasta_converter import fasta_format_converter
from  seed_finder import get_top_sequences, seed_meme_analysis, read_html_pwm, calculate_consensus
import os

def read_sequences_from_fasta_file(file_path, col_index):
    with open(file_path, 'r') as file:
        next(file)
        sequences = []
        scores = []
        for line in file:
            if line.startswith(">"):
                continue
            else:   
                if col_index == 0:
                    sequences.append(line.strip())
                else:
                    sequences.append(line.strip())
    return sequences 

def read_sequences_from_file(file_path, col_index):
    with open(file_path, 'r') as file:
        next(file)
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

def align_sequences_parallel(sequences, reference_sequence):
    substrings = generate_sorted_substrings(reference_sequence)
    tasks = [(seq, reference_sequence, substrings) for seq in sequences]
    
    with Pool(processes=cpu_count()) as pool:
        dfs = pool.map(align_sequence, tasks)

    # Combine all partial DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

def write_sequences_to_fasta(sequences, output_file):
    with open(output_file, 'w') as file:
        for i, (sequence, score) in enumerate(sequences, 1):
            file.write(f">sequence_{i}\n{sequence}\n")
       

if __name__ == "__main__":
    basepath = os.path.dirname(os.path.abspath(__file__))
    print(basepath)
    TF_list = ['Ehf','Elf2','Elf3','Elf4','Elf5','Elk1','Elk3','Elk4','Erg','Ets1','Etv1','Etv3','Etv4','Etv5','Etv6','Fli1','Gabpa','Spdef','Spic']
    for TF in TF_list[15:]:
        for col_index in [0]:
            try: 
                input_file_path = os.path.join(os.path.dirname(os.getcwd()),f'TF Binding Factors\{TF}\{TF}_8mers_top_enrichment.txt') #data file path: Uniprobe data
                sequences, scores = read_sequences_from_file(input_file_path, col_index)
                
                top_sequences = get_top_sequences(sequences, scores, top_n=20)        
                
                
                output_dir = f'MonomerBinding/{TF}/col_{col_index}'
                os.makedirs(output_dir, exist_ok=True)
                output_fasta_file = os.path.join(output_dir, f'{TF}_top_20_sequences.fasta')  #Edit this path
                
                output_file_path = os.path.join(output_dir, f'{TF}_consensus.txt')  #Edit this path
                output_csv_path = os.path.join(output_dir, f'{TF}_consensus.csv')  #Edit this path
                output_fasta_path = os.path.join(output_dir, f'{TF}_consensus.fasta')  #Edit this path

                
                write_sequences_to_fasta(top_sequences, output_fasta_file)
                
                seed_meme_analysis(output_fasta_file, col_index, 20)
                html_file = f'MonomerBinding/{TF}/col_{col_index}/Meme_of_top_20_Seeds/meme.html'  #directory path from docker
                pwm_section = read_html_pwm(html_file)
                reference_sequence = calculate_consensus(pwm_section)
                
                length_of_sequence = len(sequences[0])

                print(reference_sequence, ' ', col_index)
                if not reference_sequence:
                    raise ValueError(f"No sequence containing '{reference_sequence}' was found in the file.")

                df_aligned = align_sequences_parallel(sequences, reference_sequence)
                
                # Write to a csv file
                df_aligned.to_csv(output_csv_path, index=False)
                print(f"Data written to {output_csv_path}")
                
                columns_to_extract = [col for col in range(length_of_sequence) if col in df_aligned.columns]
                extracted_data = df_aligned[columns_to_extract]
                formatted_data = extracted_data.apply(lambda row: ''.join(row.astype(str)), axis=1)

                # Write to a plain text file
                with open(output_file_path, 'w') as file:
                    file.write(reference_sequence + '\n')
                    for line in formatted_data:
                        file.write(line + '\n')

                print(f"Data written to {output_file_path}")
                        
                # sequence_logo_generator(TF,output_file_path,col_index)
                fasta_format_converter(output_file_path, output_fasta_path)
                meme_analysis(output_fasta_path, length_of_sequence)
            
            except Exception as e:
                print(f"Error processing {TF} with col_index {col_index}: {e}")
                continue

