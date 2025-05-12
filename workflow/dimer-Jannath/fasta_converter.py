def fasta_format_converter(input_file_path, fasta_file_path):
    try:
        # Read sequences from the input file
        with open(input_file_path, 'r') as input_file:
            sequences = set()
            for line in input_file:
                sequence = line.strip()  
                if sequence: 
                    sequence = sequence.replace('-', 'N')
                    sequences.add(sequence)
        
        # Write sequences to the output FASTA file
        with open(fasta_file_path, 'w') as output_file:
            for count, sequence in enumerate(sequences, start=1):
                output_file.write(f">Sequence{count}\n")
                output_file.write(f"{sequence}\n")
        
        print(f"FASTA file created successfully at: {fasta_file_path}")
    except Exception as e:
        print(f"Error: {e}")
