import pandas as pd
from weblogo import *

def sequence_logo_generator(TF, file_path, col_index):
    fin = open(file_path)
    seqs = read_seq_data(fin)
    logodata = LogoData.from_seqs(seqs)
    logooptions = LogoOptions()
    logooptions.title = "A Logo Title"
    logooptions.color_scheme = std_color_schemes['classic']
    logooptions.yaxis_scale = 1.0
    logoformat = LogoFormat(logodata, logooptions)
    png_bytes = png_print_formatter(logodata, logoformat)

    with open(f'Working-with-TF/DimerBinding/New_SeqLogos/seed2_{TF}_{col_index}.png', 'wb') as fout:  #Edit this path
        fout.write(png_bytes)
        
        
