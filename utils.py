import re
from collections import Counter

import numpy as np
import pandas as pd
import streamlit as st


def map_locus_to_sequence_from_fasta(fasta_lines):
    locus_to_sequence_map = {}
    locus = None
    for line in fasta_lines:
        if line == "":
            continue
        elif line[0] == ">":  # new protein
            locus = line.rstrip().split(" ")[0].replace(">", "")
            description = " ".join(line.rstrip().split(" ")[1:])
            locus_to_sequence_map[locus] = {'sequence': "", 'description': description}
        else:  # protein sequence
            locus_to_sequence_map[locus]['sequence'] += line.rstrip()
    return locus_to_sequence_map


def fasta_from_locus_to_sequence_map(locus_to_sequence_map):
    lines = []
    for locus in locus_to_sequence_map:
        lines.append(f'>{locus} {locus_to_sequence_map[locus]["description"]}\n')
        lines.append(f'{locus_to_sequence_map[locus]["sequence"]}\n')
    return ''.join(lines)


def get_lines_from_uploaded_file(file):
    return file.getvalue().decode("utf-8").split("\n")


def fasta_url_download_link(tax_id, is_reviewed):
    return f'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28taxonomy_id%3A{tax_id}' \
           f'%29%20AND%20%28reviewed%3A{str(is_reviewed).lower()}%29'


def digest_sequence(sequence, enzyme_residues, min_len, max_len):
    if enzyme_residues:
        aa_to_indexes = {aa: [m.start() for m in re.finditer(aa, sequence)] for aa in enzyme_residues}
        aa_indexes = [(index, aa) for aa in aa_to_indexes for index in aa_to_indexes[aa]]
        aas = [aa for i, aa in sorted(aa_indexes, key=lambda pair: pair[0])]
        sub_sequences = re.split("|".join(enzyme_residues), sequence)

        peptides = []
        for i in range(len(sub_sequences)):
            if i < len(aas):
                peptides.append(sub_sequences[i] + aas[i])
            else:
                peptides.append(sub_sequences[i])
    else:
        peptides = []
        for i in range(len(sequence)):
            for j in range(i + 1, len(sequence)):
                if j - i < min_len:
                    continue
                if j - i > max_len:
                    break
                peptides.append(sequence[i:j])

    peptides = [peptide for peptide in peptides if min_len <= len(peptide) <= max_len]
    return peptides


def calculate_amino_acid_frequencies(proteins):
    amino_acid_counter = Counter("".join(proteins))
    amino_acids = list(amino_acid_counter.keys())
    amino_acid_counts = np.array(list(amino_acid_counter.values()))
    amino_acid_frequencies = amino_acid_counts / sum(amino_acid_counts)
    return pd.DataFrame({'Amino Acid': amino_acids, 'Frequency': amino_acid_frequencies})


def score_fasta(proteins, peptides):
    number_of_proteins = len(proteins)
    protein_average_length = np.mean([len(protein) for protein in proteins])
    number_of_peptides = len(peptides)
    number_of_unique_peptides = len(set(peptides))
    peptide_average_length = np.mean([len(peptide) for peptide in peptides])
    return number_of_proteins, protein_average_length, number_of_peptides, \
        number_of_unique_peptides, peptide_average_length


def display_stats(fasta_stats):
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Num Proteins", fasta_stats[0])
    c2.metric("Ave Protein Length", round(float(fasta_stats[1]), 2))
    c3.metric("Num Peptides", fasta_stats[2])
    c4.metric("Num Unique Peptides", fasta_stats[3])
    c5.metric("Ave Peptide Length", round(float(fasta_stats[4]), 2))
