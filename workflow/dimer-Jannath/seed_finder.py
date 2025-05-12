import os
import subprocess
from bs4 import BeautifulSoup
import json
import re

def get_top_sequences(sequences, scores, top_n=10):
    combined = list(zip(sequences, scores))
    sorted_combined = sorted(combined, key=lambda x: x[1], reverse=True)
    top_sequences = sorted_combined[:top_n]
    return top_sequences

def write_sequences_to_fasta(sequences, output_file):
    with open(output_file, 'w') as file:
        for i, (sequence, score) in enumerate(sequences, 1):
            file.write(f">sequence_{i}\n{sequence}\n")
            
def seed_meme_analysis(input_file, col_index, top_sequence_number):
        print(f"Processing file: {input_file}")
        input_abs_path = os.path.abspath(input_file)
        input_dir, input_basename = os.path.split(input_abs_path)
        print(input_basename)
        input_basename = input_basename.replace('.txt', '')
        if not os.path.exists(input_dir):
            print(f"Directory does not exist: {input_dir}")
        
        print(input_dir)
        print(f"Executing Docker command for: {input_file}")
        try:
            subprocess.run(["docker", "run", "--rm", 
                "-v", f"{input_dir}:/data",  
                "memesuite/memesuite:latest", 
                "meme", f"/data/{input_basename}", 
                "-dna", 
                "-o",
                "-nostatus",
                "-maxw", "10", 
                "-minw", "8", 
                "-nmotifs", "1", 
                "-mod", "zoops", 
                "-objfun", "classic", 
                "-revcomp", 
                "-markov_order", "0", 
                "-o", f"/data/Meme_of_top_{top_sequence_number}_Seeds"],
               check=True)
        except subprocess.CalledProcessError as e:
            print(f"Command failed with error: {e}")
            
def calculate_consensus(pwm):
    """
    Calculate the consensus sequence from a PWM.
    """
    consensus = ''
    residues = ['A', 'C', 'G', 'T']  # Assumes PWM columns are ordered as A, C, G, T
    
    for position in pwm:
        max_index = position.index(max(position))
        consensus += residues[max_index]
    
    return consensus

def read_html_pwm(file_path):
    """
    Extract PWM sections from the given list of HTML files.
    """
    with open(file_path, 'r') as file:
        html_content = file.read()
    soup = BeautifulSoup(html_content, 'html.parser')
    script_tags = soup.find_all("script")
    
    data_script = None
    for script in script_tags:
        if 'var data' in script.text:
            data_script = script.text
            break

    if data_script:
        json_str_match = re.search(r'var data = ({.*?});', data_script, re.DOTALL)
        if json_str_match:
            json_str = json_str_match.group(1)
            data = json.loads(json_str)
            if 'motifs' in data and data['motifs']:
                pwm_section = data['motifs'][0].get('pwm', [])
                if not pwm_section:
                    print(f"No PWM data found in motifs of {file_path}.")
            else:
                print(f"No 'motifs' data found in {file_path}.")
        else:
            print(f"JSON data not found in the script tag of {file_path}.")
    else:
        print(f"Script tag containing 'var data' was not found in {file_path}.")

    return pwm_section


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
    
    for j in range(seq_length):
        col_sum = sum(matrix[:, j])
        if col_sum > 0:
            matrix[:, j] /= col_sum
    
    
    information_content = np.zeros(seq_length)
    for j in range(seq_length):
        for i in range(4):
            p = matrix[i, j]
            if p > 0:
                information_content[j] -= p * np.log2(p)
    
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
    if len(motif_regions) < 2:
        # Not enough motifs found
        return None, None
    motif_regions.sort(key=len, reverse=True)
    motif_regions = motif_regions[:2]
    motif_regions.sort(key=lambda x: x[0])
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
