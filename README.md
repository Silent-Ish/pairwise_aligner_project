# Pairwise Sequence Aligner (Needleman-Wunsch)

A Python tool for performing **Global Sequence Alignment** using the Needleman-Wunsch algorithm. This tool calculates the optimal alignment between two biological sequences (DNA, RNA, or Protein) and displays the scoring matrix result alongside a visual representation of the alignment.

##  Features

* **Global Alignment:** Uses the Needleman-Wunsch dynamic programming algorithm to align sequences from end-to-end.
* **Smart Input Handling:**
    * Accepts **Raw Strings** (e.g., typing `ATCG` directly).
    * Accepts **File Paths** (FASTA files or plain text files).
    * **Drag-and-Drop Support:** Automatically cleans file paths (removes quotes and escape characters) if you drag a file into your terminal.
* **Biopython Integration:** Robust parsing of FASTA files using `Bio.SeqIO`.
* **Visual Output:** Prints alignments in readable 80-character chunks with connection lines showing matches (`|`), mismatches (`.`), and gaps (` `).

##  Prerequisites

You need Python 3 installed. This script relies on **NumPy** for matrix operations and **Biopython** for file parsing.

### Installation

1. **Clone the respiratory:**
git clone [https://github.com/Silent-Ish/pairwise_aligner_project.git](https://github.com/Silent-Ish/pairwise_aligner_project.git)
cd pairwise_aligner_project

2. **Instal dependecies**:
pip install -r requirements.txt (Alternatively: pip install numpy biopython)

```bash
pip install numpy biopython
```

##  Usage

Assuming you have named your script `aligner.py`, there are two ways to use it.

### 1. Interactive Mode (Recommended)

Simply run the script without arguments. It will prompt you for inputs. This is best if you want to drag-and-drop files into the terminal window.

```bash
python aligner.py
```

**Example Interaction:**

```text
Tip: You can drag and drop your file into the terminal window.
Please enter Sequence A (or path to file): /Users/name/data/gene1.fasta
Please enter Sequence B (or path to file): ATCCGTAGC
```

### 2. Command Line Arguments

You can pass sequences or file paths directly using flags.

```bash
# Using raw strings
python aligner.py -a ATCGTACG -b ATCGGACG

# Using file paths
python aligner.py -a data/seq1.fasta -b data/seq2.fasta
```

##  Configuration (Scoring)

The scoring system is defined at the top of the script. You can open the `.py` file and modify these constants to change how the alignment is calculated:

```python
# --- SCORING SYSTEM ---
MATCH_SCORE = 1      # Reward for a match
MISMATCH_SCORE = -1  # Penalty for a mismatch
GAP_PENALTY = -2     # Penalty for opening/extending a gap
```

##  How it Works

The tool utilizes the **Needleman-Wunsch algorithm**, a dynamic programming method.

1.  **Matrix Initialization:** Creates a grid where the axes are the two sequences.
2.  **Scoring:** Fills the grid by calculating the score for matching, mismatching, or inserting a gap based on neighbors.
3.  **Traceback:** Starts from the bottom-right corner of the matrix and follows the highest scores back to the top-left to reconstruct the optimal path.

## Author 
ismael