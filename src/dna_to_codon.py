def get_orfs(dna_sequence, min_length_aa=50):
    found_orfs = []
    seq_len = len(dna_sequence)
    stop_codons = ["TAA", "TAG", "TGA"]
    used_stops = set()
    for i in range(seq_len - 2):
        current_codon = dna_sequence[i:i+3]
        if current_codon == "ATG":
            for j in range(i, seq_len - 2, 3):
                reading_codon = dna_sequence[j:j+3]
                if reading_codon in stop_codons:
                    end_pos = j + 3
                    if end_pos not in used_stops:
                        orf_sequence = dna_sequence[i:end_pos]
                        if len(orf_sequence) >= (min_length_aa * 3):
                            found_orfs.append({
                                "start_position": i + 1,
                                "end_position": end_pos,
                                "sequence": orf_sequence
                            })
                            used_stops.add(end_pos)
                    break
    return found_orfs