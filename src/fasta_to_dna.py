import re
from typing import Dict

def fasta_to_dna(fasta_file: str) -> Dict[str, str]:
    """
    Parses a FASTA file and extracts the DNA sequences.
    
    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        Dict[str, str]: A dictionary mapping sequence IDs to their corresponding DNA sequences.
    
    Raises:
        ValueError: If a sequence contains non-IUPAC characters.
    """
    sequences = {}
    current_id = "Unnamed_Sequence"
    iupac_pattern = re.compile(r'[^ATGCN]') 
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current_id = line[1:].split()[0]
                sequences[current_id] = []
            else:
                if current_id not in sequences:
                    sequences[current_id] = []
                cleaned_line = line.upper()

                if iupac_pattern.search(cleaned_line):
                    invalid_chars = set(iupac_pattern.findall(cleaned_line))
                    raise ValueError(f"Sequence '{current_id}' contains invalid non-IUPAC characters: {invalid_chars}")
                sequences[current_id].append(cleaned_line)
                
    for seq_id in sequences:
        sequences[seq_id] = "".join(sequences[seq_id])
    return sequences
