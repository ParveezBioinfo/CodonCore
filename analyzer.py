import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

def get_basic_stats(sequence_string):
    """Calculates sequence length and overall GC content."""
    seq = Seq(sequence_string)
    length = len(seq)
    gc_percent = gc_fraction(seq) * 100
    return length, round(gc_percent, 2)

def calculate_sliding_window_gc(sequence_string, window_size=50):
    """Calculates GC content in small chunks to plot the line graph."""
    positions = []
    gc_values = []
    
    # Slide across the sequence
    for i in range(0, len(sequence_string) - window_size + 1, 10):
        window_seq = sequence_string[i:i+window_size]
        gc = gc_fraction(Seq(window_seq)) * 100
        positions.append(i)
        gc_values.append(round(gc, 2))
        
    return positions, gc_values

def find_orfs(sequence_string, min_length):
    """Finds Open Reading Frames (ORFs) in all 6 frames."""
    orfs = []
    seq = Seq(sequence_string)
    
    # 3 forward frames and 3 reverse frames
    frames = [
        (1, seq), (2, seq[1:]), (3, seq[2:]),
        (-1, seq.reverse_complement()), 
        (-2, seq.reverse_complement()[1:]), 
        (-3, seq.reverse_complement()[2:])
    ]
    
    for frame_id, frame_seq in frames:
        protein = frame_seq.translate(to_stop=False)
        peptides = str(protein).split("*")
        
        current_pos = 0
        for peptide in peptides:
            start_idx = peptide.find("M") # Methionine is the start codon
            
            if start_idx != -1:
                orf_seq = peptide[start_idx:]
                orf_length = len(orf_seq) * 3
                
                if orf_length >= min_length:
                    orfs.append({
                        "Frame": frame_id,
                        "Start Pos": (current_pos + start_idx) * 3 + 1,
                        "Stop Pos": (current_pos + start_idx) * 3 + orf_length,
                        "Length (bp)": orf_length,
                        "Translated Amino Acid Sequence": orf_seq
                    })
            current_pos += len(peptide) + 1 # Move past the stop codon
            
    return pd.DataFrame(orfs)