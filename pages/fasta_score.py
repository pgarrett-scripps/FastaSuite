import re
from collections import Counter

import numpy as np
import pandas as pd
import streamlit as st
from matplotlib import pyplot as plt
from scipy import spatial

from decoy_util import VALID_AMINO_ACIDS
from utils import map_locus_to_sequence_from_fasta

st.title("FASTA Score")
with st.expander("Help"):
    st.markdown("""
        Scores a FASTA file according to decoy/target statistics.
        The "ideal" FASTA file contains target and decoy proteins which statistically resemble each other.
        
        Enzyme Sites: The digestion sites. Use K,R for trypsin. 
        
        Min/Max Peptide Lengths: Keep only digested peptides within these bounds (bounds are inclusive)
        
        Decoy Flag: The flag used to represent Decoy Proteins
    """)

fasta_file = st.file_uploader("Choose a fasta file", type=".fasta")
enzyme_residues = st.multiselect("Enzyme Sites", list(VALID_AMINO_ACIDS), ['K', 'R'], help="Residues to cleave after. For Trypsin use (K & R)")
min_len, max_len = st.slider("Min/Max Peptide Lengths", 0, 100, [6, 50], help="Sets range of supported peptides (min <= len_of_peptide <= max)")
decoy_flag = st.text_input("Decoy Flag", "DECOY_", help="Flag used to identify Decoy peptides (if fasta contains '>DECOY_sp|XXXX|YYYY' set to 'DECOY_')")

if st.button("Run"):
    if not fasta_file:
        st.warning('Upload a FASTA file!')

    if fasta_file is not None:
        fasta_lines = fasta_file.getvalue().decode("utf-8").split("\n")
        locus_to_sequence_map = map_locus_to_sequence_from_fasta(fasta_lines)

        target_sequences, decoy_sequences = [], []
        for locus in locus_to_sequence_map:
            if decoy_flag in locus:
                decoy_sequences.append(locus_to_sequence_map[locus]['sequence'])
            else:
                target_sequences.append(locus_to_sequence_map[locus]['sequence'])

        def digest_sequence(sequence, enzyme_residues, min_len, max_len):
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

            peptides = [peptide for peptide in peptides if min_len <= len(peptide) <= max_len]
            return peptides

        target_peptides = [peptide for sequence in target_sequences for peptide in digest_sequence(sequence, enzyme_residues, min_len, max_len)]
        decoy_peptides = [peptide for sequence in decoy_sequences for peptide in digest_sequence(sequence, enzyme_residues, min_len, max_len)]


        def display_fasta_stats(target_proteins, target_peptides, decoy_proteins, decoy_peptides):
            st.subheader("Target")

            c1, c2, c3, c4, c5 = st.columns(5)
            num_target_proteins = len(target_proteins)
            ave_target_protein_len = round(np.mean([len(seq) for seq in target_proteins]),2)
            num_target_peptides = len(target_peptides)
            num_unique_target_peptides = len(set(target_peptides))
            ave_target_peptide_length = round(np.mean([len(seq) for seq in target_peptides]),2)

            target_amino_acid_frequencies = Counter("".join(target_proteins))
            aas = list(target_amino_acid_frequencies.keys())
            counts = np.array(list(target_amino_acid_frequencies.values()))
            freqs = counts/sum(counts)
            target_df = pd.DataFrame({'Amino Acid':aas, 'Frequency': freqs})
            target_frequncy_score = round(1 - spatial.distance.cosine(target_df.Frequency, target_df.Frequency),4)

            c1.metric("Num Proteins", num_target_proteins)
            c2.metric("Ave Protein Length", ave_target_protein_len)
            c3.metric("Num Peptides", num_target_peptides)
            c4.metric("Num Unique Peptides", num_unique_target_peptides)
            c5.metric("Ave Peptide Length", ave_target_peptide_length)


            st.subheader("Decoy")
            c1, c2, c3, c4, c5 = st.columns(5)
            num_decoy_proteins = len(decoy_proteins)
            ave_decoy_protein_len = round(np.mean([len(seq) for seq in decoy_proteins]),2)
            num_decoy_peptides = len(decoy_peptides)
            num_decoy_target_peptides = len(set(decoy_peptides))
            ave_decoy_peptide_length = round(np.mean([len(seq) for seq in decoy_peptides]),2)

            decoy_amino_acid_frequencies = Counter("".join(decoy_proteins))
            aas = list(decoy_amino_acid_frequencies.keys())
            counts = np.array(list(decoy_amino_acid_frequencies.values()))
            freqs = counts/sum(counts)
            decoy_df = pd.DataFrame({'Amino Acid':aas, 'Frequency': freqs})

            df = pd.merge(target_df, decoy_df, on='Amino Acid', how='outer', suffixes=(' Target', ' Decoy'))
            decoy_frequncy_score = round(1 - spatial.distance.cosine(df['Frequency Target'], df['Frequency Decoy']),4)


            c1.metric("Num Proteins", num_decoy_proteins, num_decoy_proteins- num_target_proteins)
            c2.metric("Ave Protein Length", ave_decoy_protein_len, ave_decoy_protein_len - ave_target_protein_len)
            c3.metric("Num Peptides", num_decoy_peptides, num_decoy_peptides-num_target_peptides)
            c4.metric("Num Unique Peptides", num_decoy_target_peptides, num_decoy_target_peptides-num_unique_target_peptides)
            c5.metric("Ave Peptide Length", ave_decoy_peptide_length, ave_decoy_peptide_length-ave_target_peptide_length)

            st.subheader("Target/Decoy")
            c1, c2, = st.columns(2)
            c1.metric("AA Frequency Score", decoy_frequncy_score, decoy_frequncy_score-target_frequncy_score)
            c2.metric("Shared Peptides", len(set(target_peptides).intersection(set(decoy_peptides))))

            with st.expander("Amino Acid Frequency"):

                fig, ax = plt.subplots()
                indexes = list(range(len(df['Amino Acid'].values)))
                ax.bar(indexes, df['Frequency Target'].values, alpha=0.5, label="Target")
                ax.bar(indexes, df['Frequency Decoy'].values, alpha =0.5, label="Decoy")
                ax.set_xlabel('Amino Acids')
                ax.set_title('Amino Acid Frequency')
                ax.set_xticks(indexes)
                ax.set_xticklabels(df['Amino Acid'].values)
                ax.legend()
                st.pyplot(fig)
                df.set_index('Amino Acid', inplace=True)

                st.table(df.round(5))

        display_fasta_stats(target_sequences, target_peptides, decoy_sequences, decoy_peptides)

