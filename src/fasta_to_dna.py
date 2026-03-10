def fasta_to_dna(fasta_file):
    sequences = {}
    current_id = "Unnamed_Sequence"
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current_id = line[1:]
                sequences[current_id] = []
            else:
                if current_id not in sequences:
                    sequences[current_id] = []
                sequences[current_id].append(line.upper())
    for seq_id in sequences:
        sequences[seq_id] = "".join(sequences[seq_id])
    return sequences