import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random

def generate_example_sequences(output_file, 
                               motif1="ATGAC", 
                               motif2="GTCAT",
                               num_sequences=100,
                               min_spacing=0,
                               max_spacing=10,
                               sequence_length=20,
                               mutation_rate=0.1):
    """
    Generate example sequences containing dimer motifs with variable spacing.
    
    Args:
        output_file: Path to save the generated sequences
        motif1: First motif sequence
        motif2: Second motif sequence
        num_sequences: Number of sequences to generate
        min_spacing: Minimum spacing between motifs
        max_spacing: Maximum spacing between motifs
        sequence_length: Total length of each sequence
        mutation_rate: Probability of introducing mutations
    """
    nucleotides = ['A', 'T', 'G', 'C']
    sequences = []
    scores = []
    
    for i in range(num_sequences):
        # Determine spacing for this sequence
        spacing = random.randint(min_spacing, max_spacing)
        
        # Potentially mutate motifs
        motif1_mut = ''.join([n if random.random() > mutation_rate else random.choice(nucleotides)
                             for n in motif1])
        motif2_mut = ''.join([n if random.random() > mutation_rate else random.choice(nucleotides)
                             for n in motif2])
        
        # Create spacing segment
        spacing_segment = ''.join([random.choice(nucleotides) for _ in range(spacing)])
        
        # Calculate how much flanking sequence we need
        total_motif_length = len(motif1_mut) + spacing + len(motif2_mut)
        remaining_length = sequence_length - total_motif_length
        
        if remaining_length < 0:
            # Handle case where motifs + spacing exceed sequence length
            # Just use motifs with reduced spacing
            spacing_segment = spacing_segment[:max(0, spacing + remaining_length)]
            sequence = motif1_mut + spacing_segment + motif2_mut
            # Trim if still too long
            sequence = sequence[:sequence_length]
        else:
            # Add flanking sequences
            left_flank_size = random.randint(0, remaining_length)
            right_flank_size = remaining_length - left_flank_size
            
            left_flank = ''.join([random.choice(nucleotides) for _ in range(left_flank_size)])
            right_flank = ''.join([random.choice(nucleotides) for _ in range(right_flank_size)])
            
            sequence = left_flank + motif1_mut + spacing_segment + motif2_mut + right_flank
        
        # Calculate a mock "binding score" - higher when motifs are less mutated
        motif1_accuracy = sum(1 for a, b in zip(motif1, motif1_mut) if a == b) / len(motif1)
        motif2_accuracy = sum(1 for a, b in zip(motif2, motif2_mut) if a == b) / len(motif2)
        score = 0.5 * (motif1_accuracy + motif2_accuracy)
        
        # Add some noise to the score
        score = score * (0.8 + 0.4 * random.random())
        score = min(1.0, score)
        
        sequences.append(sequence)
        scores.append(score)
    
    # Write sequences to file
    with open(output_file, 'w') as f:
        f.write("Sequence\tReverseComplement\tScore\n")
        for seq, score in zip(sequences, scores):
            # Generate reverse complement
            rev_comp = ''.join({'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[n] for n in reversed(seq))
            f.write(f"{seq}\t{rev_comp}\t{score:.4f}\n")
    
    print(f"Generated {num_sequences} example sequences with dimer motifs.")
    print(f"Saved to {output_file}")
    
    # Also write a FASTA version
    fasta_file = output_file.replace('.txt', '.fasta')
    write_sequences_to_fasta(zip(sequences, scores), fasta_file)
    
    return sequences, scores

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
                file.write(f">sequence_{i}_score_{score:.4f}\n{sequence}\n")
            else:
                file.write(f">sequence_{i}\n{item}\n")
    
    print(f"Sequences written to FASTA file: {output_file}")

if __name__ == "__main__":
    # Create a directory for our example data
    os.makedirs("example_data", exist_ok=True)
    
    # Generate example sequences with dimer motifs and variable spacing
    output_file = "example_data/example_dimer_sequences.txt"
    
    # Define the motifs
    motif1 = "ATGAC"
    motif2 = "GTCAT"
    
    # Generate sequences with the dimer motifs and variable spacing
    sequences, scores = generate_example_sequences(
        output_file,
        motif1=motif1,
        motif2=motif2,
        num_sequences=100,
        min_spacing=0,
        max_spacing=10,
        sequence_length=20,
        mutation_rate=0.1
    )
    
    print(f"Generated example sequences with dimer motifs {motif1} and {motif2}")
    print(f"Variable spacing between 0 and 10 nucleotides")
    print(f"Example files created: {output_file} and {output_file.replace('.txt', '.fasta')}")