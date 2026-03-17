import argparse
from pathlib import Path
from dna_to_rna import dna_to_rna
from fasta_to_dna import fasta_to_dna
from dna_to_codon import get_orfs
from dna_to_protein import translate_sequence
from results_export import export_orfs_to_txt
from sequence_properties import calculate_dna_properties, calculate_rna_properties, calculate_protein_properties

base_path = Path(__file__).resolve().parent.parent

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SeqProfiler: DNA analysis tool.")
    parser.add_argument("--min-length", type=int, default=None, help="Minimum ORF size (in amino acids) [default: 50]")
    parser.add_argument("--input", type=str, default=None, help="Path to input directory containing .fasta files [default: data/]")
    parser.add_argument("--output", type=str, default=None, help="Path to output directory for results [default: results/]")
    args = parser.parse_args()

    print("=== SeqProfiler (BATCH MODE) ===\n")
    
    if args.min_length is None:
        min_length = 50
        print(f"No minimum size specified. Using default: {min_length}")
    else:
        min_length = args.min_length
        print(f"Minimum ORF size: {min_length}")

    data_dir = Path(args.input).resolve() if args.input else base_path / "data"
    results_dir = Path(args.output).resolve() if args.output else base_path / "results"
    
    print(f"Input directory: {data_dir}")
    print(f"Output directory: {results_dir}")
    print("-" * 44)

    if not data_dir.exists():
        print(f"Error: Input directory '{data_dir}' does not exist.")
        exit(1)
    
    results_dir.mkdir(parents=True, exist_ok=True)
    
    fasta_files = list(data_dir.glob("*.fasta"))
    
    if not fasta_files:
        print("No .fasta files found in the data/ directory.")
    else:
        print(f"{len(fasta_files)} file(s) found. Starting analysis...\n")
        for file_path in fasta_files:
            print(f"Processing {file_path.name}...")
            sequences_dict = fasta_to_dna(file_path)
            all_file_orfs = []
            
            for seq_id, main_sequence in sequences_dict.items():
                seq_len = len(main_sequence)
                reverse_complement = main_sequence.translate(str.maketrans("ATCG", "TAGC"))[::-1]
                positive_orfs = get_orfs(main_sequence, min_length_aa=min_length)
                for orf in positive_orfs:
                    orf["sequence_id"] = seq_id
                    orf["strand"] = "Forward"
                    p3, p1 = translate_sequence(orf["sequence"])
                    orf["protein_3l"] = p3
                    orf["protein_1l"] = p1
                    orf["rna"] = dna_to_rna(orf["sequence"])
                    orf["dna_props"] = calculate_dna_properties(orf["sequence"])
                    orf["rna_props"] = calculate_rna_properties(orf["rna"])
                    orf["prot_props"] = calculate_protein_properties(p1)
                    
                negative_orfs = get_orfs(reverse_complement, min_length_aa=min_length)
                for orf in negative_orfs:
                    orf["sequence_id"] = seq_id
                    orf["strand"] = "Reverse"
                    p3, p1 = translate_sequence(orf["sequence"])
                    orf["protein_3l"] = p3
                    orf["protein_1l"] = p1
                    true_start = seq_len - orf["end_position"] + 1
                    true_end = seq_len - orf["start_position"] + 1
                    orf["start_position"] = true_start
                    orf["end_position"] = true_end
                    orf["rna"] = dna_to_rna(orf["sequence"])
                    orf["dna_props"] = calculate_dna_properties(orf["sequence"])
                    orf["rna_props"] = calculate_rna_properties(orf["rna"])
                    orf["prot_props"] = calculate_protein_properties(p1)
                    
                all_file_orfs.extend(positive_orfs + negative_orfs)
                
            all_file_orfs = sorted(all_file_orfs, key=lambda x: (x["sequence_id"], x["start_position"]))
            output_name = f"results_{file_path.stem}.txt"
            output_path = results_dir / output_name
            export_orfs_to_txt(all_file_orfs, output_path)
            
        print("\nAll analyses completed successfully.")