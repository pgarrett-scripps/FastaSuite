import pandas as pd
import streamlit as st
from matplotlib import pyplot as plt
from scipy import spatial

from constants import VALID_AMINO_ACIDS, FASTA_SCORE_HELP_MESSAGE, DECOY_FLAG_HELP_MESSAGE, \
    MIN_MAX_PEPTIDE_LENGTH_HELP_MESSAGE, ENZYME_HELP_MESSAGE
from utils import map_locus_to_sequence_from_fasta, digest_sequence, score_fasta, calculate_amino_acid_frequencies, \
    display_stats

st.title("Score FASTA")
with st.expander("Help"):
    st.markdown(FASTA_SCORE_HELP_MESSAGE)

fasta_file = st.file_uploader(label="Choose a fasta file", type=".fasta")
enzyme_residues = st.multiselect(label="Enzyme Sites", options=list(VALID_AMINO_ACIDS), default=['K', 'R'],
                                 help=ENZYME_HELP_MESSAGE)
min_len, max_len = st.slider(label="Min/Max Peptide Lengths", min_value=0, max_value=100, value=[6, 50],
                             help=MIN_MAX_PEPTIDE_LENGTH_HELP_MESSAGE)
decoy_flag = st.text_input(label="Decoy Flag", value="DECOY_", help=DECOY_FLAG_HELP_MESSAGE)

MAX_PEPTIDES = 250_000
with st.expander('Advance Params'):
    max_peptides = st.number_input(label='maximum number of peptides to digest', value=MAX_PEPTIDES)

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

        target_peptides = []
        for sequence in target_sequences:
            peptides = digest_sequence(sequence, enzyme_residues, min_len, max_len)
            target_peptides.extend(peptides)
            if len(target_peptides) >= max_peptides:
                target_peptides = target_peptides[:max_peptides]
                break

        decoy_peptides = []
        for sequence in decoy_sequences:
            peptides = digest_sequence(sequence, enzyme_residues, min_len, max_len)
            decoy_peptides.extend(peptides)
            if len(decoy_peptides) >= max_peptides:
                decoy_peptides = decoy_peptides[:max_peptides]
                break

        st.subheader("Target")
        fasta_target_stats = score_fasta(target_sequences, target_peptides)
        display_stats(fasta_target_stats)

        st.subheader("Decoy")
        fasta_decoy_stats = score_fasta(decoy_sequences, decoy_peptides)
        display_stats(fasta_decoy_stats)

        target_amino_acid_frequency_df = calculate_amino_acid_frequencies(target_sequences)
        decoy_amino_acid_frequency_df = calculate_amino_acid_frequencies(decoy_sequences)

        amino_acid_frequency_df = pd.merge(target_amino_acid_frequency_df, decoy_amino_acid_frequency_df,
                                           on='Amino Acid', how='outer', suffixes=(' Target', ' Decoy'))
        amino_acid_frequency_score = 1 - spatial.distance.cosine(amino_acid_frequency_df['Frequency Target'],
                                                                 amino_acid_frequency_df['Frequency Decoy'])

        st.subheader("Target/Decoy")
        c1, c2, = st.columns(2)
        c1.metric("AA Frequency Score", round(amino_acid_frequency_score, 4),
                  round(amino_acid_frequency_score - 1, 4))
        c2.metric("Shared Peptides", len(set(target_peptides).intersection(set(decoy_peptides))))

        with st.expander("Amino Acid Frequency"):
            fig, ax = plt.subplots()
            indexes = list(range(len(amino_acid_frequency_df['Amino Acid'].values)))
            ax.bar(indexes, amino_acid_frequency_df['Frequency Target'].values, alpha=0.5, label="Target")
            ax.bar(indexes, amino_acid_frequency_df['Frequency Decoy'].values, alpha=0.5, label="Decoy")
            ax.set_xlabel('Amino Acids')
            ax.set_title('Amino Acid Frequency')
            ax.set_xticks(indexes)
            ax.set_xticklabels(amino_acid_frequency_df['Amino Acid'].values)
            ax.legend()
            st.pyplot(fig)
            amino_acid_frequency_df.set_index('Amino Acid', inplace=True)
            st.table(amino_acid_frequency_df.round(5))
