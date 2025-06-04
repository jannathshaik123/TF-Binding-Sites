# TFBS Anchor Residue Identification Project

A comprehensive bioinformatics toolkit for identifying and analyzing anchor residues within Transcription Factor Binding Sites (TFBS) to advance targeted drug design.

## 📋 Table of Contents

- [Overview](#-overview)
- [Features](#-features)
- [Installation](#-installation)
- [File Structure](#-file-structure)
- [Usage](#-usage)
- [BINGO Algorithm](#-bingo-algorithm)
- [Requirements](#--requirements)
- [Examples](#-examples)
- [Contributing](#-contributing)

## 🎯 Overview

This project aims to identify critical anchor residues within transcription factor binding sites by analyzing UniProbe data using computational methods[7]. The ultimate goal is to enable the development of novel therapeutic strategies that can inhibit gene expression in diseased states by targeting these critical residues.

The project implements three different methodologies, with the **BINGO Algorithm** being the most promising approach for accurate anchor residue identification.

## ✨ Features

- **Multi-Method Analysis**: Three distinct approaches for anchor residue identification
- **MEME Integration**: Leverages MEME Suite tools through Docker containers
- **Sequence Alignment**: Advanced alignment algorithms for dataset comparison
- **Motif Discovery**: Automated motif identification and analysis
- **Visualization**: Sequence logo generation and PWM analysis
- **Format Conversion**: Seamless conversion between different sequence formats

**Note**: This workflow implementation is adapted from the original repository by [Secret-Ambush](https://github.com/Secret-Ambush/Identifying-Anchor-Residues-in-TFBS).

## 🚀 Installation

### Prerequisites

- Python 3.9.13 or higher
- Docker (version 25.0.3 or higher)
- Required Python packages (see requirements section)

### Setup Steps

1. **Clone the repository**

```bash
git clone https://github.com/jannathshaik123/TF-Binding-Sites
cd TF-Binding-Sites-main/workflow
```

2. **Install Python dependencies**

```bash
pip install pandas numpy tqdm beautifulsoup4 weblogo
```

3. **Pull MEME Suite Docker image**

```bash
docker pull memesuite/memesuite:latest
```

4. **Prepare your data**
   - Place your UniProbe `*_8mers.txt` files in the appropriate directory
   - Ensure proper directory structure for TF data

## 📁 File Structure

```
project/
├── alignment_utils.py      # Core alignment functions and substring generation
├── fasta_converter.py      # Format conversion utilities
├── main.py                # Main processing pipeline and entry point
├── seed_finder.py         # Top sequence identification and MEME analysis
├── seq_logo_gen.py        # Sequence logo generation utilities
├── using_meme.py          # MEME Suite integration functions
└── README.md              # This file
```

### Core Files Description

**`alignment_utils.py`**[1]

- Contains the `generate_sorted_substrings()` function for creating ordered substring lists
- Implements the `align_sequence()` function for sequence alignment with reference sequences
- Handles DataFrame operations for storing alignment results

**`fasta_converter.py`**[2]

- Provides `fasta_format_converter()` function to convert text files to FASTA format
- Replaces gaps (-) with unknown nucleotides (N) for compatibility
- Essential for MEME Suite input preparation

**`main.py`**[3]

- Central processing hub with command-line interface
- Implements the complete pipeline from data reading to result generation
- Supports both single file and directory batch processing
- Handles parallel sequence alignment for improved performance

**`seed_finder.py`**[4]

- Extracts top sequences based on enrichment scores
- Integrates with MEME Suite through Docker containers
- Parses HTML output to extract Position Weight Matrices (PWM)
- Calculates consensus sequences from PWM data

**`seq_logo_gen.py`**[5]

- Generates visual sequence logos using WebLogo
- Creates publication-ready graphics for motif visualization

**`using_meme.py`**[6]

- Handles MEME Suite Docker container execution
- Configures motif discovery parameters
- Manages output file organization

## 🔧 Usage

### Basic Usage

**Process a single TF file:**

```bash
python main.py --file path/to/GATA4_8mers.txt --col 0 --top 20
```

**Process multiple files in a directory:**

```bash
python main.py --dir path/to/tf_directory --top 20
```

**Process specific TF list:**

```bash
python main.py --tf_list GATA4 ETV5 --col 0
```

### Command Line Arguments

- `--file`: Path to a single sequence file
- `--dir`: Directory containing multiple sequence files
- `--col`: Column index for sequence data (0 or 1, default: 0)
- `--top`: Number of top sequences to analyze (default: 20)
- `--pattern`: File pattern to match when scanning directories
- `--tf_list`: List of TF names to process

## 🎯 BINGO Algorithm

The most promising approach that eliminates positioning uncertainty through comprehensive dataset alignment.

### Introduction

The BINGO (Binding Site Genomic Organization) Algorithm represents the most advanced method in this toolkit. Unlike the previous methods, it doesn't use binning and instead aligns the entire dataset against the top-scoring motif[7].

### How BINGO Works

#### Step 1: Seed Selection

```python
# Takes the highest scoring 8-mer as the reference seed
seed = top_sequence_from_dataset  # e.g., "AGATAAGG"
```

#### Step 2: Substring Generation

The algorithm generates ALL possible substrings from the seed sequence, ordered by length (longest first) and alphabetically:

```python
def generate_sorted_substrings(seed):
    # For seed "AGATAAGG", generates:
    # 8-mer: AGATAAGG
    # 7-mer: AGATAAGG, GATAAGG
    # 6-mer: AGATAA, GATAAGG, ATAAGG
    # ... and so on
    substrings = [seed[i:j] for i in range(len(seed))
                  for j in range(i + 1, len(seed) + 1)]
    return sorted(substrings, key=lambda x: (-len(x), x))
```

#### Step 3: Intelligent Sequence Alignment

For each sequence in the dataset, BINGO:

1. **Searches for the longest matching substring** with the seed
2. **Calculates alignment offset** based on substring positions
3. **Aligns sequences** in a standardized coordinate system
4. **Handles mismatches gracefully** by finding the best possible alignment

```python
def align_sequence(sequence, seed, substrings):
    # Find longest common substring
    for substring in substrings:
        if substring in sequence:
            # Calculate alignment offset
            seq_pos = sequence.find(substring)
            seed_pos = seed.find(substring)
            offset = seed_pos - seq_pos

            # Align and return positioned sequence
            return align_with_offset(sequence, offset)
```

#### Step 4: Anchor Residue Identification

After alignment, BINGO:

- Creates a position-specific frequency matrix
- Identifies highly conserved positions across all sequences
- Determines anchor residues based on conservation patterns
- Generates visual sequence logos for interpretation

### Key Advantages of BINGO

1. **Position Accuracy**: Eliminates ambiguity in anchor residue positioning
2. **Comprehensive Analysis**: Uses entire dataset rather than binned subsets
3. **Flexible Alignment**: Handles sequence variations intelligently
4. **Visual Output**: Generates clear sequence logos for result interpretation

### BINGO Algorithm Implementation

The algorithm is implemented across two main functions in `alignment_utils.py`:

**Algorithm 1: Sequence Reading and Substring Generation**

```python
def read_sequences_from_file(file_path):
    # Extracts sequences, removes gaps (.)
    # Returns clean sequence list

def generate_sorted_substrings(seed):
    # Creates all possible substrings
    # Sorts by length (desc) then alphabetically
    # Returns unique substring list
```

**Algorithm 2: Sequence Alignment and Data Extraction**

```python
def align_sequences(sequences, seed):
    # Generates substrings from seed
    # Creates alignment DataFrame
    # For each sequence:
    #   - Finds best matching substring
    #   - Calculates alignment offset
    #   - Places sequence in correct position
    # Returns aligned DataFrame
```

### Expected Output

BINGO produces:

- **Aligned sequence file** (`.txt` format)
- **FASTA formatted results** for further analysis
- **Sequence logos** showing conservation patterns
- **Position-specific anchor residue identification**

## 💻 Requirements

### System Requirements

- **OS**: macOS Sonoma 14.5 (tested), Linux, Windows with WSL
- **Memory**: 16 GB RAM recommended
- **Processor**: Multi-core CPU for parallel processing

### Python Dependencies

```
pandas>=1.3.0
numpy>=1.21.0
tqdm>=4.62.0
beautifulsoup4>=4.10.0
weblogo>=3.7.0
multiprocessing (built-in)
subprocess (built-in)
argparse (built-in)
```

### External Tools

- **Docker**: For MEME Suite containerized execution
- **MEME Suite**: Motif discovery and analysis (via Docker)

## 📊 Examples

### Example 1: GATA4 Analysis

```bash
python main.py --file GATA4_8mers_top_enrichment.txt --col 0 --top 20
```

**Expected Output:**

```
Processing monomer GATA4 from GATA4_8mers_top_enrichment.txt
Reference sequence for GATA4: AGATAA
Data written to MonomerBinding/GATA4/col_0/GATA4_consensus.csv
FASTA file created successfully at: MonomerBinding/GATA4/col_0/GATA4_consensus.fasta
```

### Example 2: Batch Processing

```bash
python main.py --dir ./TF_data/ --top 50
```

### Example 3: Multiple TF Analysis

```bash
python main.py --tf_list GATA4 ETV5 NR3C1 --col 0
```

## 🧪 Validation

The project has been tested on well-characterized transcription factors:

- **GATA4**: Key cardiac development regulator[7]
- **ETV5**: Important for stem cell maintenance and differentiation[7]

Results can be validated against:

- Crystal structure data from PDB
- Published anchor residue studies
- Experimental binding affinity measurements

## 🤝 Contributing

We welcome contributions to improve the toolkit:

1. **Fork the repository**
2. **Create a feature branch** (`git checkout -b feature/new-method`)
3. **Commit changes** (`git commit -am 'Add new analysis method'`)
4. **Push to branch** (`git push origin feature/new-method`)
5. **Create Pull Request**

### Areas for Contribution

- Additional motif discovery algorithms
- Enhanced visualization tools
- Support for additional file formats
- Performance optimizations
- Machine learning integration for anchor prediction

## 📝 Citation

If you use this toolkit in your research, please cite:

```
Goswami, R., Jothi, J.A.A., & Ghoshdastidar, D. (2025).
Identification and Analysis of Anchor Residues within Transcription Factor Binding Sites (TFBS).
Computational Biology Research Project.
```

For questions, issues, or feature requests:

- **Open an issue** on GitHub

---

**Note**: This toolkit is designed for research purposes in computational biology and drug design. Results should be validated through experimental methods before therapeutic applications.
