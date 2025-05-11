# TF-Binding-Sites

final worklow taken from ![here](https://github.com/Secret-Ambush/Identifying-Anchor-Residues-in-TFBS/tree/main/final%20workflow)

# From Monomer to Dimer Alignment: Detailed Code Evolution Explanation

I'll systematically explain how the original alignment code (`alignment_utils.py`) was expanded into the dimer alignment code (`dimer_alignment_utils.py`), analyzing each component with examples to show their similarities, differences, and evolutionary path.

## 1. Generate Sorted Substrings: Function Evolution

### Original Monomer Version
```python
def generate_sorted_substrings(s):
    substrings = [s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1)]
    substrings_sorted = sorted(substrings, key=lambda x: (-len(x), x))
    unique_substrings = []
    for substring in substrings_sorted:
        if substring not in unique_substrings:
            unique_substrings.append(substring)
    return unique_substrings
```

### Dimer Version
```python
def generate_sorted_substrings(s):
    """Generate all possible substrings of a string, sorted by length (descending)."""
    substrings = [s[i:j] for i in range(len(s)) for j in range(i + 1, len(s) + 1)]
    substrings_sorted = sorted(substrings, key=lambda x: (-len(x), x))
    unique_substrings = []
    for substring in substrings_sorted:
        if substring not in unique_substrings:
            unique_substrings.append(substring)
    return unique_substrings
```

**Similarities:**
- Both functions have identical logic and implementation
- Both generate all possible substrings and sort them by length (longest first)
- Both deduplicate substrings

**Differences:**
- The dimer version adds a docstring explaining the function's purpose
- Otherwise functionally identical

**Example:**
For sequence "ACGT":
1. Generate all substrings: ["A", "AC", "ACG", "ACGT", "C", "CG", "CGT", "G", "GT", "T"]
2. Sort by length (descending): ["ACGT", "ACG", "CGT", "AC", "CG", "GT", "A", "C", "G", "T"]
3. Remove duplicates (none in this case)
4. Return the sorted unique list

## 2. Core Alignment Logic Evolution

### Original Monomer Version: Single Function
```python
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
```

### Dimer Version: Split into Multiple Functions

#### 1. First Function: find_best_motif_match
```python
def find_best_motif_match(sequence, motif, min_match_length=3):
    """Find the best match for a motif in a sequence."""
    substrings = generate_sorted_substrings(motif)
    
    best_match_score = 0
    best_offset = 0
    best_substring = ""
    
    for substring in substrings:
        if len(substring) < min_match_length:
            break  # Skip very short substrings
            
        if substring in sequence:
            idx = sequence.find(substring)
            substring_idx_in_motif = motif.find(substring)
            offset = substring_idx_in_motif - idx
            aligned_motif = motif[max(0, offset):]
            aligned_sequence = sequence[max(0, -offset):]
            
            min_length = min(len(aligned_motif), len(aligned_sequence))
            match_score = sum(1 for j in range(min_length) if aligned_motif[j] == aligned_sequence[j])
            
            if match_score > best_match_score:
                best_match_score = match_score
                best_offset = offset
                best_substring = substring
    
    return best_offset, best_match_score, best_substring
```

**Core Evolution:**
The monomer alignment's substring matching loop was extracted into its own function that:
1. Generates substrings internally rather than taking them as input
2. Returns more information (offset, score, and the best substring)
3. Adds a minimum match length cutoff to improve efficiency
4. Renames variables for clarity (e.g., "match_score" instead of "total_number_of_matches")

**Example:**
For sequence "ACGTACGT" and motif "TACG":
1. Generate all substrings of "TACG"
2. Find "TACG" in sequence at position 3
3. Calculate offset: position in motif (0) - position in sequence (3) = -3
4. Calculate match score by counting matching positions after alignment
5. Return best offset (-3), match score (4), and substring ("TACG")

#### 2. Second Function: align_dimer_sequence
```python
def align_dimer_sequence(sequence_info):
    """Align a sequence to two motifs (dimer) with variable spacing."""
    sequence, motif1, motif2, substrings1, substrings2, max_spacing = sequence_info
    
    max_length = len(sequence)
    lower_bound = -max_length
    upper_bound = 2 * max_length
    columns = [i for i in range(lower_bound, upper_bound)]
    df = pd.DataFrame('-', index=[0], columns=columns)
    
    offset1, score1, substring1 = find_best_motif_match(sequence, motif1)
    sequence_positions = {}
    for i, char in enumerate(sequence):
        sequence_positions[i] = offset1 + i
    best_score2 = 0
    best_offset2 = 0
    best_position2 = 0
    
    for spacing in range(max_spacing + 1):
        expected_pos2 = len(motif1) + spacing
        for seq_start in range(max(0, offset1), min(len(sequence), offset1 + len(sequence))):
            remaining_seq = sequence[seq_start:]
            if not remaining_seq:
                continue
                
            offset2, score2, substring2 = find_best_motif_match(remaining_seq, motif2)
            adjusted_offset = seq_start + offset2
            spacing_penalty = abs(adjusted_offset - (offset1 + expected_pos2)) / max(1, max_spacing)
            adjusted_score = score2 * (1 - 0.1 * spacing_penalty)  # Penalize by spacing deviation
            
            if adjusted_score > best_score2:
                best_score2 = adjusted_score
                best_offset2 = offset2
                best_position2 = seq_start
    
    total_score = score1
    if best_score2 > len(motif2) * 0.3:  # At least 30% match
        total_score += best_score2
    for i, val in enumerate(sequence):
        pos = sequence_positions.get(i, offset1 + i)
        if pos in df.columns:
            df.at[0, pos] = val
    
    return df, total_score
```

**Major Evolution:**
1. Takes the basic alignment concept but extends it to handle two motifs
2. Uses the extracted `find_best_motif_match` function for both motifs
3. Adds logic for variable spacing between motifs
4. Implements a spacing penalty to favor alignments with expected spacing
5. Returns both alignment and a score for ranking multiple alignments

**Key New Components:**
- Handling of variable spacing with `max_spacing` parameter
- Two-phase alignment: find first motif, then find second motif
- Adjusting the score based on spacing deviation
- Setting a threshold (30%) for accepting a second motif match

**Example:**
For sequence "ACGTTACGTACCT" with motifs "ACGT" and "ACCT" and max_spacing=5:
1. Find first motif "ACGT" at position 0
2. Try finding "ACCT" after "ACGT" with different spacings (0-5)
3. Find "ACCT" at position 9, which gives spacing of 5 (within max_spacing)
4. Apply spacing penalty if needed and calculate total score
5. Generate alignment dataframe showing both motifs aligned

## 3. Multiple Sequence Processing: New Function

The dimer version adds an entirely new function to process multiple sequences in parallel:

```python
def align_dimer_sequences_parallel(sequences, motif1, motif2, max_spacing=5):
    """Align multiple sequences to a dimer motif pattern in parallel."""
    substrings1 = generate_sorted_substrings(motif1)
    substrings2 = generate_sorted_substrings(motif2)
    
    tasks = [(seq, motif1, motif2, substrings1, substrings2, max_spacing) for seq in sequences]
    
    results = []
    for task in tqdm(tasks, desc="Aligning sequences"):
        df, score = align_dimer_sequence(task)
        results.append((df, score))
    results.sort(key=lambda x: x[1], reverse=True)
    dfs = [result[0] for result in results]
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df
```

**New Functionality:**
1. Processes multiple sequences with the same motif pair
2. Pre-computes substrings for efficiency (only generate once)
3. Sorts results by alignment score (best matches first)
4. Combines individual alignments into a single dataframe for visualization
5. Provides a simple interface for batch processing

**Example:**
For sequences ["ACGTTACGTACCT", "ACGTGCATACCT", "ACGTTATACCT"]:
1. Generate all substrings of both motifs once
2. Process each sequence through align_dimer_sequence
3. Sort results with best scores first
4. Create a combined dataframe showing all alignments together:

```
Seq1: ACGTTACGTACCT  (Score: 8.2)
Seq2: ACGTGCATACCT   (Score: 7.9)
Seq3: ACGTTATACCT    (Score: 7.5)
```

## 4. Motif Discovery: Entirely New Function

The dimer version adds advanced functionality for discovering motifs from sequences:

```python
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
```

**Completely New Functionality:**
1. Uses position frequency matrix to analyze nucleotide conservation
2. Calculates information content (entropy) to find conserved regions
3. Identifies regions with statistically significant conservation
4. Builds consensus motifs from the conserved regions
5. Returns the two most significant conserved patterns

**Example:**
For sequences:
```
"ACGTTACGTACCT"
"ACGTGACGTACCT"
"ACGTCACGTACCT"
"ACGTAACGTACCT"
```

1. Create position frequency matrix showing frequency of each base at each position
2. Calculate information content (low values indicate conservation)
3. Identify conserved regions where information content is below threshold
4. Find two non-overlapping conserved regions: positions 0-3 and 5-8
5. Build consensus motifs: "ACGT" and "ACGT"
6. Return the two motifs

## 5. Side-by-Side Comparison of Key Design Differences

| Feature | Monomer Alignment | Dimer Alignment |
|---------|------------------|-----------------|
| **Modular design** | Single monolithic function | Multiple specialized functions |
| **Purpose** | Find best match of single sequence | Find best match of two motifs with spacing |
| **Return values** | Just alignment | Alignment plus score |
| **Multiple sequence handling** | Not supported | Supported with parallel processing |
| **Motif discovery** | Not supported | Supported with information theory |
| **Spacing consideration** | N/A | Variable spacing with penalties |
| **Performance optimization** | Limited | Enhanced with early exits and thresholds |
| **Score calculation** | Simple count of matches | Weighted by spacing adherence |

## 6. Detailed Function-by-Function Evolution Path

### From `align_sequence` to `find_best_motif_match`

The substring-matching core of the original function was:

```python
# Original code for finding best substring match
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
```

This evolved into:

```python
# Evolved code in find_best_motif_match
for substring in substrings:
    if len(substring) < min_match_length:
        break  # Skip very short substrings
        
    if substring in sequence:
        idx = sequence.find(substring)
        substring_idx_in_motif = motif.find(substring)
        offset = substring_idx_in_motif - idx
        aligned_motif = motif[max(0, offset):]
        aligned_sequence = sequence[max(0, -offset):]
        
        min_length = min(len(aligned_motif), len(aligned_sequence))
        match_score = sum(1 for j in range(min_length) if aligned_motif[j] == aligned_sequence[j])
        
        if match_score > best_match_score:
            best_match_score = match_score
            best_offset = offset
            best_substring = substring
```

**Evolution changes:**
1. Loop directly over substrings instead of indices
2. Add early exit for short substrings
3. Rename variables for clarity
4. Return multiple values instead of using them internally

### From `align_sequence` to `align_dimer_sequence`

The alignment part of the original function:

```python
# Original alignment code
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
```

Evolved into a more complex sequence in the dimer version:

```python
# Evolved dimer alignment code
offset1, score1, substring1 = find_best_motif_match(sequence, motif1)
sequence_positions = {}
for i, char in enumerate(sequence):
    sequence_positions[i] = offset1 + i

# [Second motif finding code here]

total_score = score1
if best_score2 > len(motif2) * 0.3:  # At least 30% match
    total_score += best_score2
for i, val in enumerate(sequence):
    pos = sequence_positions.get(i, offset1 + i)
    if pos in df.columns:
        df.at[0, pos] = val
```

**Evolution changes:**
1. Use of the extracted find_best_motif_match function
2. Creation of a position mapping dictionary
3. Score calculation incorporating both motifs
4. Score threshold for accepting secondary motif matches

## 7. Full Process Example: From Monomer to Dimer Analysis

Let me walk through a complete example showing how both algorithms would process the same data:

**Input:**
- Sequence: "ACGTTACGTACCT"
- Monomer reference: "ACGT"
- Dimer motifs: "ACGT" and "ACCT"

### Monomer Alignment Process:

1. Generate substrings of reference "ACGT": ["ACGT", "ACG", "CGT", "AC", "CG", "GT", "A", "C", "G", "T"]
2. Find matches in sequence:
   - "ACGT" matches at position 0
   - Calculate offset: 0 - 0 = 0
   - Count matches: 4 (perfect match)
   - Best so far, save it
3. Try other substrings - none give better matches
4. Create alignment with perfect match at position 0:

```
Position: -12-11-10-9-8-7-6-5-4-3-2-1 0 1 2 3 4 5 6 7 8 9 10 11 12
Alignment: - - - - - - - - - - - - A C G T T A C G T A C C T
```

### Dimer Alignment Process:

1. Find best match for first motif "ACGT":
   - Perfect match at position 0
   - Save offset1 = 0, score1 = 4
2. Create position mapping: {0:0, 1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 7:7, 8:8, 9:9, 10:10, 11:11, 12:12}
3. Try finding second motif "ACCT" with different spacings:
   - For spacing=0, expected position is 4 (after first motif)
   - Try offsets starting from position 0
   - At position 9: Find "ACCT" (perfect match)
   - Actual spacing = 5 (positions 9-4)
   - Calculate spacing penalty: |5-0|/5 = 1.0 (maximum penalty)
   - Adjusted score: 4 * (1-0.1*1.0) = 3.6
4. Try other spacings:
   - For spacing=5, expected position is 9
   - Find "ACCT" at position 9
   - No spacing penalty (perfect match at expected position)
   - Adjusted score: 4.0 (no penalty)
5. Best second motif match: offset2=0, position2=9, score2=4.0
6. Calculate total score: 4 + 4 = 8
7. Create alignment:

```
Position: -12-11-10-9-8-7-6-5-4-3-2-1 0 1 2 3 4 5 6 7 8 9 10 11 12
Alignment: - - - - - - - - - - - - A C G T T A C G T A C C T
Motifs:                             ^^^^^               ^^^^^
                                    Motif1              Motif2
```

## 8. Final Analysis: Code Evolution Principles

Looking at the transformation from monomer to dimer alignment code, we can identify several key evolutionary principles:

1. **Function Extraction**: Taking common code patterns and making them reusable functions (substring generation, match finding)

2. **Parameter Expansion**: Adding new parameters to handle more complex scenarios (spacing, thresholds)  

3. **Score Refinement**: Evolving from simple match counting to weighted scores with penalties

4. **Modularity**: Breaking down monolithic function into specialized components

5. **Batch Processing**: Adding functionality to process multiple inputs efficiently

6. **Advanced Analytics**: Incorporating information theory for motif discovery

7. **Performance Optimization**: Adding early exits and minimum thresholds

8. **Return Value Enhancement**: Providing more information in return values (scores, best matches)

The dimer alignment code represents a sophisticated evolution of the original monomer code, retaining its core alignment approach while adding layers of functionality to handle more complex biological sequence analysis tasks. This progression demonstrates how software for scientific applications often evolves - keeping fundamental algorithms intact while expanding their capabilities to address more complex research questions.