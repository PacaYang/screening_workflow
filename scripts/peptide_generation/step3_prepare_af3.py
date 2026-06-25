#!/usr/bin/env python3
"""Prepare AF3 input JSON files for protein-peptide complexes.

The target may be a multi-chain oligomer (e.g. the TL1A homotrimer in 2O0O).
Each target chain is emitted as a separate AF3 protein chain, and the designed
peptide is appended as the final chain.
"""

import json
import os
import argparse
import string
from pathlib import Path


def parse_fasta(fasta_file):
    """Parse FASTA file and return list of (header, sequence) tuples (order preserved)."""
    sequences = []
    current_header = None
    current_seq = []

    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))

    return sequences


def create_af3_json(target_chains, peptide_seq, name, output_path):
    """Create AF3 JSON for a multi-chain target plus a designed peptide.

    target_chains: list of target chain sequences (one entry per chain).
    peptide_seq:   the designed peptide sequence (appended as the last chain).

    The AF3 jobs run with --norun_data_pipeline, so every protein chain must
    carry explicit (empty) MSA and template fields. Empty unpaired/paired MSAs
    put AF3 into single-sequence mode for that chain.
    """
    chain_ids = list(string.ascii_uppercase)
    sequences = []

    def protein_chain(chain_id, seq):
        return {
            "protein": {
                "id": chain_id,
                "sequence": seq,
                "unpairedMsa": "",
                "pairedMsa": "",
                "templates": [],
            }
        }

    for i, seq in enumerate(target_chains):
        sequences.append(protein_chain(chain_ids[i], seq))

    # Peptide gets the next available chain letter (e.g. D for a trimer target)
    peptide_id = chain_ids[len(target_chains)]
    sequences.append(protein_chain(peptide_id, peptide_seq))

    data = {
        "name": name,
        "dialect": "alphafold3",
        "version": 1,
        "sequences": sequences,
        "modelSeeds": [99],
    }

    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)


def main():
    parser = argparse.ArgumentParser(description="Generate AF3 input JSONs for peptide validation")
    parser.add_argument('--protein-fasta', required=True,
                        help="Target protein FASTA file (monomer sequence)")
    parser.add_argument('--peptide-dir', required=True,
                        help="Directory with ProteinMPNN output FASTA files")
    parser.add_argument('--output-dir', required=True, help="Output directory for JSON files")
    parser.add_argument('--protein-name', required=True, help="Protein name for output")
    parser.add_argument('--target-copies', type=int, default=1,
                        help="Number of target chains (e.g. 3 for a homotrimer). "
                             "The monomer sequence is replicated this many times.")
    parser.add_argument('--no-skip-native', action='store_true', default=False,
                        help="Include the first FASTA entry instead of skipping it. "
                             "Use when input is not ProteinMPNN output (e.g. BoltzGen sequences).")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Parse target monomer sequence. If the FASTA already contains multiple
    # records, treat each as a distinct chain; otherwise replicate the single
    # monomer --target-copies times.
    protein_seqs = parse_fasta(args.protein_fasta)
    if len(protein_seqs) == 1:
        monomer = protein_seqs[0][1]
        target_chains = [monomer] * args.target_copies
    else:
        target_chains = [seq for _, seq in protein_seqs]

    print(f"Target: {len(target_chains)} chain(s), "
          f"each {len(target_chains[0])} aa")

    # Process all peptide FASTA files
    peptide_dir = Path(args.peptide_dir)
    fasta_files = sorted(peptide_dir.glob('*.fa'))

    json_count = 0
    for fasta_file in fasta_files:
        entries = parse_fasta(fasta_file)
        backbone_id = fasta_file.stem

        # ProteinMPNN output: first entry is the native input sequence, the
        # remaining entries are the designed samples. Skip the native one,
        # unless --no-skip-native is set (e.g. for BoltzGen sequences).
        design_entries = entries if args.no_skip_native else entries[1:]

        for idx, (header, seq) in enumerate(design_entries, 1):
            # ProteinMPNN concatenates chains with '/'. The backbone has the
            # target monomer and the designed peptide; the peptide is the last
            # chain segment.
            chains = seq.split('/')
            peptide_seq = chains[-1].strip()

            json_name = f"{args.protein_name}_{backbone_id}_seq{idx:03d}"
            json_path = Path(args.output_dir) / f"{json_name}.json"
            create_af3_json(target_chains, peptide_seq, json_name, json_path)
            json_count += 1

    print(f"Created {json_count} AF3 input JSON files in {args.output_dir}")


if __name__ == "__main__":
    main()
