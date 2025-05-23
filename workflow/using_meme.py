import os
import subprocess

def meme_analysis(input_file, max_seq_length):
        print(f"Processing file: {input_file} ***")
        input_abs_path = os.path.abspath(input_file)
        print(f"Processing file: {input_abs_path} ***")
        input_dir, input_basename = os.path.split(input_abs_path)
        print(f"Base name: {input_basename} ***")
        if not os.path.exists(input_dir):
            print(f"Directory does not exist: {input_dir}")
        
        output_basename = input_basename.replace('.fasta', '')

        print(f"Executing Docker command for: {input_file} ***")
        try:
            subprocess.run(["docker", "run", "--rm", 
                "-v", f"{input_dir}:/data",  
                "memesuite/memesuite:latest", 
                "meme", f"/data/{input_basename}", 
                "-dna", 
                "-o",
                "-nostatus",
                "-maxw", f"{max_seq_length}", 
                "-minw", "8", 
                "-nmotifs", "1", 
                "-mod", "zoops", 
                "-objfun", "classic", 
                "-revcomp", 
                "-markov_order", "0", 
                "-o", f"/data/{output_basename}"],
               check=True)
        except subprocess.CalledProcessError as e:
            print(f"Command failed with error: {e}")
