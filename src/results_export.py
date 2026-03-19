from typing import List, Dict

def export_orfs_to_txt(orfs: List[Dict], output_file: str) -> None:
    """
    Exports ORF analysis results to a structured text file.
    
    Args:
        orfs (List[Dict]): List of ORF data dictionaries.
        output_file (str): Path to the output file.
    """
    with open(output_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("SeqProfiler: ORF Analysis Report\n")
        f.write("=" * 60 + "\n\n")
        
        if not orfs:
            f.write("No ORFs found matching the criteria.\n")
            return

        for orf in orfs:
            f.write(f"> Sequence ID: {orf['sequence_id']}\n")
            f.write(f"Strand: {orf['strand']}\n")
            f.write(f"Position: {orf['start_position']} - {orf['end_position']} (bp)\n")
            f.write("-" * 40 + "\n")
            
            # DNA Properties
            dp = orf['dna_props']
            f.write(f"[DNA]  Length: {dp['length']} bp | Mass: {dp['mass_da']:.2f} Da | GC: {dp['gc_prop']:.1f}% | Tm: {dp['tm']:.1f} °C\n")
            
            # RNA Properties
            rp = orf['rna_props']
            f.write(f"[RNA]  Mass: {rp['mass_da']:.2f} Da\n")
            
            # Protein Properties & Sequence
            pp = orf['prot_props']
            f.write(f"[PROT] Length: {len(orf['protein_1l'])} aa | Mass: {pp['mass_kda']:.2f} kDa | pI: {pp['pi']:.2f} | Ext. Coeff: {pp['ext_coeff']}\n")
            f.write(f"\n[PROTEIN SEQUENCE (1L)]\n{orf['protein_1l']}\n")
            f.write("\n" + "." * 60 + "\n\n")
