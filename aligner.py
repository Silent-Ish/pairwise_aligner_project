import numpy as np
import argparse
import os
from Bio import SeqIO

# --- SCORING SYSTEM ---
# You can change these to see how alignment changes!
MATCH_SCORE = 1
MISMATCH_SCORE = -1
GAP_PENALTY = -2

def read_sequence(input_argument):
    """
    Smart reader:
    1. Cleans input (removes quotes from drag-and-drop).
    2. Handles terminal escape characters (like backslashes before spaces).
    3. Checks if the input is a file path.
    4. If it is a file, reads it (supports FASTA or raw text).
    5. If it's not a file, treats the input as the raw sequence string itself.
    """
    if not input_argument:
        return ""

    # Remove quotes that might appear if you drag-and-drop a file path
    clean_path = input_argument.strip('"').strip("'")
    
    # FIX: Remove the backslash that terminals add before spaces (e.g. "Science\ /Bio...")
    clean_path = clean_path.replace(r'\ ', ' ')

    if os.path.exists(clean_path):
        try:
            # Try parsing as FASTA first (Best Practice)
            record = SeqIO.read(clean_path, "fasta")
            return str(record.seq).upper()
        except:
            # If not FASTA, just read the raw text from the file
            with open(clean_path, 'r') as f:
                # Remove newlines and whitespace
                return f.read().replace('\n', '').strip().upper()
    else:
        # It's not a file, so it must be the raw sequence (e.g., "ATCG")
        return input_argument.upper()

def needleman_wunsch(seq1, seq2):
    """
    Implementation of the Needleman-Wunsch Global Alignment Algorithm.
    """
    n = len(seq1)
    m = len(seq2)
    
    # 1. Initialize the Scoring Matrix (Size: n+1 x m+1)
    # We use NumPy because it's faster than Python lists
    score_matrix = np.zeros((n + 1, m + 1))
    
    # 2. Fill the first row and column with Gap Penalties
    # (Because aligning sequence to "nothing" costs points)
    for i in range(n + 1):
        score_matrix[i][0] = i * GAP_PENALTY
    for j in range(m + 1):
        score_matrix[0][j] = j * GAP_PENALTY
        
    # 3. Fill the rest of the matrix
    # We look at 3 neighbors: Diagonal (Match/Mismatch), Top (Gap), Left (Gap)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            
            # Calculate scores
            if seq1[i-1] == seq2[j-1]:
                match = score_matrix[i-1][j-1] + MATCH_SCORE
            else:
                match = score_matrix[i-1][j-1] + MISMATCH_SCORE
            
            delete = score_matrix[i-1][j] + GAP_PENALTY
            insert = score_matrix[i][j-1] + GAP_PENALTY
            
            # Take the best score
            score_matrix[i][j] = max(match, delete, insert)

    # 4. Traceback (Reconstruct the alignment path)
    align1 = ""
    align2 = ""
    i, j = n, m
    
    while i > 0 and j > 0:
        score_current = score_matrix[i][j]
        score_diagonal = score_matrix[i-1][j-1]
        score_up = score_matrix[i][j-1]
        score_left = score_matrix[i-1][j]
        
        # Check if we came from Diagonal (Match/Mismatch)
        if seq1[i-1] == seq2[j-1]:
            is_match = MATCH_SCORE
        else:
            is_match = MISMATCH_SCORE

        if score_current == score_diagonal + is_match:
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif score_current == score_left + GAP_PENALTY:
            align1 += seq1[i-1]
            align2 += "-"
            i -= 1
        elif score_current == score_up + GAP_PENALTY:
            align1 += "-"
            align2 += seq2[j-1]
            j -= 1

    # Finish remaining gaps if we hit the top or left wall
    while i > 0:
        align1 += seq1[i-1]
        align2 += "-"
        i -= 1
    while j > 0:
        align1 += "-"
        align2 += seq2[j-1]
        j -= 1

    return align1[::-1], align2[::-1], score_matrix[n][m]

def main():
    parser = argparse.ArgumentParser(description="Pairwise Sequence Aligner (Needleman-Wunsch)")
    # UPDATED: Arguments are now optional (required=False)
    parser.add_argument("-a", "--seq1", help="First Sequence (File or String)", required=False)
    parser.add_argument("-b", "--seq2", help="Second Sequence (File or String)", required=False)
    
    args = parser.parse_args()
    
    # NEW LOGIC: If arguments are missing, ask for them interactively
    if not args.seq1:
        print("Tip: You can drag and drop your file into the terminal window.")
        args.seq1 = input("Please enter Sequence A (or path to file): ").strip()
        
    if not args.seq2:
        args.seq2 = input("Please enter Sequence B (or path to file): ").strip()

    # Read the inputs using our new helper function
    s1 = read_sequence(args.seq1)
    s2 = read_sequence(args.seq2)
    
    print(f"\nAligning Sequence A ({len(s1)} bp) and Sequence B ({len(s2)} bp)...")
    
    # Run the alignment
    a_aligned, b_aligned, score = needleman_wunsch(s1, s2)
    
    # Print the results nicely
    print("\n--- Alignment Result ---")
    
    # Create the connection lines (| for match, . for mismatch)
    connections = ""
    for k in range(len(a_aligned)):
        if a_aligned[k] == b_aligned[k]:
            connections += "|"
        elif a_aligned[k] == "-" or b_aligned[k] == "-":
            connections += " "
        else:
            connections += "."
            
    # Print in chunks of 80 characters so it fits on screen
    chunk_size = 80
    for k in range(0, len(a_aligned), chunk_size):
        a_chunk = a_aligned[k:k+chunk_size]
        b_chunk = b_aligned[k:k+chunk_size]
        conn_chunk = connections[k:k+chunk_size]
        
        print(f"SeqA: {a_chunk}")
        print(f"      {conn_chunk}")
        print(f"SeqB: {b_chunk}")
        print("")

    print(f"Final Score: {score}")
    print("------------------------\n")

if __name__ == "__main__":
    main()