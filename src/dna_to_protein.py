genetic_code = {
    "AAA": "Lys", "AAT": "Asn", "AAG": "Lys", "AAC": "Asn",
    "ATA": "Ile", "ATT": "Ile", "ATG": "Met", "ATC": "Ile",
    "AGA": "Arg", "AGT": "Ser", "AGG": "Arg", "AGC": "Ser",
    "ACA": "Thr", "ACT": "Thr", "ACG": "Thr", "ACC": "Thr",
    "TAA": "Ter", "TAT": "Tyr", "TAG": "Ter", "TAC": "Tyr",
    "TTA": "Leu", "TTT": "Phe", "TTG": "Leu", "TTC": "Phe",
    "TGA": "Ter", "TGT": "Cys", "TGG": "Trp", "TGC": "Cys",
    "TCA": "Ser", "TCT": "Ser", "TCG": "Ser", "TCC": "Ser",
    "GAA": "Glu", "GAT": "Asp", "GAG": "Glu", "GAC": "Asp",
    "GTA": "Val", "GTT": "Val", "GTG": "Val", "GTC": "Val",
    "GGA": "Gly", "GGT": "Gly", "GGG": "Gly", "GGC": "Gly",
    "GCA": "Ala", "GCT": "Ala", "GCG": "Ala", "GCC": "Ala",
    "CAA": "Gln", "CAT": "His", "CAG": "Gln", "CAC": "His",
    "CTA": "Leu", "CTT": "Leu", "CTG": "Leu", "CTC": "Leu",
    "CGA": "Arg", "CGT": "Arg", "CGG": "Arg", "CGC": "Arg",
    "CCA": "Pro", "CCT": "Pro", "CCG": "Pro", "CCC": "Pro"
}

def dna_to_protein(dna_codon):
    return genetic_code.get(dna_codon, "X")

aa_3_to_1 = {
    "Lys": "K", "Asn": "N", "Ile": "I", "Met": "M", "Arg": "R", "Ser": "S",
    "Thr": "T", "Ter": "*", "Tyr": "Y", "Leu": "L", "Phe": "F", "Cys": "C",
    "Trp": "W", "Glu": "E", "Asp": "D", "Val": "V", "Gly": "G", "Ala": "A",
    "Gln": "Q", "His": "H", "Pro": "P", "X": "X"
}

def translate_sequence(dna_sequence):
    protein_3l = ""
    protein_1l = ""
    for i in range(0, len(dna_sequence) - 3, 3):
        codon = dna_sequence[i:i+3]
        aa_3 = dna_to_protein(codon)
        if aa_3 != "Ter":
            protein_3l += aa_3
            protein_1l += aa_3_to_1.get(aa_3, "?") 
    return protein_3l, protein_1l