## Overview

**ORF_Detector** is a purely foundational, from-scratch bioinformatics tool written in Python. It parses standard FASTA files, identifies valid Open Reading Frames (ORFs), transcribes DNA into mRNA, translates sequences into amino acids, and calculates key biochemical properties (Mass, Tm, pI, Extinction Coefficient) for all identified molecules.

## Project Status

**Completed / Pedagogical.** This repository was created as a foundational learning exercise and portfolio piece. It is fully functional for its intended purpose, but it may or may not receive active updates or maintenance in the future.

## The "Vanilla Python" Philosophy

If you are reviewing this code, you might wonder why I manually hardcoded a 64-codon dictionary, built custom string-slicing loops, or wrote binary search algorithms to calculate isoelectric points instead of simply importing `Bio.Seq` from Biopython.

**This is an intentional and pedagogical constraint.** As a 3rd-year biology student in Grenoble, my goal was to deeply understand and computationally replicate the mechanical logic of the ribosomal machinery and the mathematics behind biochemical properties. Relying on "black box" libraries would defeat the purpose of the exercise. Every algorithmic choice was written explicitly to prove a granular understanding of genomic translation and molecular biophysics.

## Core Features

- **Custom Multi-FASTA Parsing:** Efficiently handles multi-line sequence concatenation and whitespace stripping without external dependencies.

- **Strict ORF Detection:** Scans both forward (+) and reverse (-) strands for the `ATG` start codon, stepping exactly in-frame (modulo 3) until a stop codon (`TAA`, `TAG`, `TGA`) is reached.

- **Intelligent Redundancy Filtering:** Utilizes a heuristic `used_stops` tracking system to prioritize the longest functional reading frame and prevent reporting smaller, nested ORFs.

- **Biochemical Profiling:** Calculates physical and chemical properties for all sequences from scratch. This includes DNA/RNA exact mass, sequence GC/AT content, melting temperature (Tm), and protein molecular weight (kDa), theoretical isoelectric point (pI), and extinction coefficient.

- **Bi-Directional Coordinate Mapping:** Accurately maps the 1-based start/end coordinates of reverse-complement ORFs back to the original forward strand reference.

- **Batch Processing:** Automatically processes all `.fasta` files in the `/data` directory and exports structured, detailed reports to `/results`.

## Architecture

The directory structure intentionally mirrors the biological flow of information:

```text
ORF_Detector/
├── data/                       # Input .fasta files
├── results/                    # Output text reports
└── src/
    ├── dna_to_codon.py         # Frame reading & Start/Stop identification
    ├── dna_to_protein.py       # Translation (1-letter & 3-letter formats)
    ├── dna_to_rna.py           # Transcription logic
    ├── fasta_to_dna.py         # Sequence extraction & multi-FASTA mapping
    ├── main.py                 # Main execution point
    ├── results_export.py       # Structured report generation
    └── sequence_properties.py  # Biochemical property calculations (Mass, pI, Tm, etc.)
```

## Usage

1. **Prepare Data:** Place your target `.fasta` files into the `data/` directory.
    
2. **Run Pipeline:** Execute the main script from the root directory:

```bash
python src/main.py
```

3. **Configure:** Enter your desired minimum ORF size in amino acids (defaults to 50) when prompted.

4. **Review:** Retrieve detailed analysis reports in the `results/` folder.


## License

This project is licensed under the MIT License.