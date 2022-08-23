import random
from pathlib import Path

import streamlit as st

from utils import map_locus_to_sequence_from_fasta, fasta_from_locus_to_sequence_map

fasta_file = st.file_uploader("Choose a fasta file", type=".fasta")

decoy_strategy = st.radio("Decoy Strategy", ('reverse', 'shuffle'), index=0)
decoy_flag = st.text_input("Decoy Flag", "DECOY_")

random_seed = None
if decoy_strategy == 'shuffle':
    random_seed = st.number_input("random number seed", value=7878)

if st.button("Generate Decoys"):

    if not fasta_file:
        st.warning('Upload a FASTA file!')

    if fasta_file is not None:

        with st.spinner("Generating decoys..."):

            target_protein_col, decoy_protein_col = st.columns(2)
            fasta_lines = fasta_file.getvalue().decode("utf-8").split("\n")

            locus_to_sequence_map = map_locus_to_sequence_from_fasta(fasta_lines)
            decoy_locus_to_sequence_map = {}
            for locus in locus_to_sequence_map:
                decoy_locus = decoy_flag + locus
                description = locus_to_sequence_map[locus]['description']

                sequence = locus_to_sequence_map[locus]['sequence']
                if decoy_strategy == 'reverse':
                    decoy_sequence = sequence[::-1]
                else:
                    decoy_sequence = ''.join(random.sample(sequence, len(sequence)))
                decoy_locus_to_sequence_map[decoy_locus] = {'description':description, 'sequence':decoy_sequence}

            target_protein_col.metric(label="Target Proteins", value=len(locus_to_sequence_map))
            decoy_protein_col.metric(label="Decoy Proteins", value=len(decoy_locus_to_sequence_map))

            new_fasta_lines = fasta_from_locus_to_sequence_map({**locus_to_sequence_map, **decoy_locus_to_sequence_map})
            st.download_button("Download FASTA", "".join(new_fasta_lines), file_name=f"{Path(fasta_file.name).stem}_filter.fasta")

            st.balloons()
