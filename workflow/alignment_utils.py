import pandas as pd
from tqdm import tqdm

def generate_sorted_substrings(s):
    substrings = [s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1)]
    substrings_sorted = sorted(substrings, key=lambda x: (-len(x), x))
    unique_substrings = []
    for substring in substrings_sorted:
        if substring not in unique_substrings:
            unique_substrings.append(substring)
    return unique_substrings

def align_sequence(sequence_info):
    sequence, reference_sequence, substrings = sequence_info
    max_length = len(sequence)
    lower_bound = 1 - max_length
    upper_bound = 2 * max_length - 1
    columns = [i for i in range(lower_bound, upper_bound)]
    df = pd.DataFrame('-', index=[0], columns=columns)

    best_substring_index = -1
    prev_best_number_of_matches = 0
    best_offset = 0

    # Iterate over substrings to find the best alignment
    for i in tqdm(range(len(substrings)), desc="Finding best substring match"):
        if substrings[i] in sequence:
            current_substring = substrings[i]
            idx = sequence.find(current_substring)
            substring_index_in_seed = reference_sequence.find(current_substring)
            offset = substring_index_in_seed - idx

            aligned_seed = reference_sequence[max(0, offset):]
            aligned_sequence = sequence[max(0, -offset):]

            min_length = min(len(aligned_seed), len(aligned_sequence))
            total_number_of_matches = sum(1 for j in range(min_length) if aligned_seed[j] == aligned_sequence[j])

            if total_number_of_matches > prev_best_number_of_matches:
                prev_best_number_of_matches = total_number_of_matches
                best_substring_index = i
                best_offset = offset

    if best_substring_index != -1:
        best_substring = substrings[best_substring_index]
        idx = sequence.find(best_substring)
        substring_index_in_seed = reference_sequence.find(best_substring)
        offset = best_offset

        # Place aligned data in the respective indices
        indices = range(offset, offset + len(sequence))
        for i, val in zip(indices, sequence):
            if i in df.columns: 
                df.at[0, i] = val
        
    return df
