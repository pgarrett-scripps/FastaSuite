from pathlib import Path

import streamlit as st

from decoy_util import VALID_AMINO_ACIDS, reverse_sequence, static_shuffle_sequence, shuffle_sequence, \
    exchange_sequence, shift_reverse_sequence, make_sequence_markov_model, make_locus_name_markov_model, \
    make_gene_name_markov_model, predict_sequence_from_markov_model, predict_locus_name_from_markov_model, \
    predict_gene_name_from_markov_model
from utils import map_locus_to_sequence_from_fasta, fasta_from_locus_to_sequence_map

with st.expander("Help"):
    st.markdown("""

    This app will generate decoy proteins sequences given a fasta file.

    **Input**

    **FASTA File**: FASTA file containing target proteins
    
    **Decoy Flag:** The flag to prepend locus names ("Reverse_" -> ">Reverse_sp|XXXX|YYYY")
    
    **Decoy Strategies**
    
    - **reverse:** reverse the sequence
    - **shuffle:** randomly shuffle the sequence
    - **static shuffle:** randomly shuffle residues between given static residues
    - **markov:** trains a markov chain model on target proteins and randomly generates decoys with length between min/max
    - **exchange:** exchange the provided residues with the replacement
    - **shifted reversal:** reverse the sequence, then switch the given amino acids with their predecessor 

    """)

fasta_file = st.file_uploader("Choose a fasta file", type=".fasta")

decoy_flag = st.text_input("Decoy Flag", "DECOY_")
decoy_strategy = st.radio("Decoy Strategy", ('reverse', 'shuffle', 'markov', 'exchange', 'shifted reversal'), index=0)

random_seed = None
if decoy_strategy == 'shuffle' or decoy_strategy == 'markov':
    random_seed = st.number_input("random number seed", value=7878)

static_amino_acids = None
if decoy_strategy == 'shuffle':
    static_amino_acids = st.multiselect(label='Static Residues', options=list(VALID_AMINO_ACIDS), default=None)

aa_exchange_map = {}
if decoy_strategy == 'exchange':
    with st.expander("Residue Map"):
        for aa in VALID_AMINO_ACIDS:
            aa_exchange_map[aa] = st.selectbox(label=f'{aa}', options=list(VALID_AMINO_ACIDS), index=list(VALID_AMINO_ACIDS).index(aa))

shifted_amino_acids = None
if decoy_strategy == 'shifted reversal':
    shifted_amino_acids = set(st.multiselect("Shifted Residues", options=list(VALID_AMINO_ACIDS),
                                             default=['K', 'R']))

if st.button("Generate Decoys"):

    if not fasta_file:
        st.warning('Upload a FASTA file!')

    if fasta_file is not None:

        with st.spinner("Generating decoys..."):

            target_protein_col, decoy_protein_col = st.columns(2)
            fasta_lines = fasta_file.getvalue().decode("utf-8").split("\n")

            locus_to_sequence_map = map_locus_to_sequence_from_fasta(fasta_lines)
            lengths = [len(locus_to_sequence_map[locus]['sequence']) for locus in locus_to_sequence_map]
            min_len, max_len = min(lengths), max(lengths)
            org_name = list(locus_to_sequence_map.keys())[0].split("|")[2].split("_")[1]

            if decoy_strategy == 'markov':
                sequences = [locus_to_sequence_map[locus]['sequence'] for locus in locus_to_sequence_map]
                sequence_model = make_sequence_markov_model(sequences)

                locus_names = [locus.split("|")[1] for locus in locus_to_sequence_map]
                locus_name_model = make_locus_name_markov_model(locus_names)

                gene_names = [locus.split("|")[1].split("_")[0] for locus in locus_to_sequence_map]
                gene_name_model = make_gene_name_markov_model(gene_names)

            decoy_locus_to_sequence_map = {}
            for locus in locus_to_sequence_map:
                decoy_locus = decoy_flag + locus
                decoy_description = locus_to_sequence_map[locus]['description']

                sequence = locus_to_sequence_map[locus]['sequence']
                if decoy_strategy == 'reverse':
                    decoy_sequence = reverse_sequence(sequence)
                elif decoy_strategy == 'shuffle':
                    if static_amino_acids:
                        decoy_sequence = static_shuffle_sequence(sequence, static_amino_acids)
                    else:
                        decoy_sequence = shuffle_sequence(sequence)
                elif decoy_strategy == 'markov':

                    decoy_sequence = predict_sequence_from_markov_model(sequence_model, min_len, max_len)
                    locus_name = predict_locus_name_from_markov_model(locus_name_model, 2, 10)
                    gene_name = predict_gene_name_from_markov_model(gene_name_model, 2, 10)

                    if locus_name is None:
                        locus_name = "XXXX"

                    if gene_name is None:
                        gene_name = "YYYY_ORG"

                    decoy_description = "Made up Protein"
                    decoy_locus = f"{decoy_flag}|{locus_name}|{gene_name}"

                elif decoy_strategy == 'exchange':
                    decoy_sequence = exchange_sequence(sequence, aa_exchange_map)

                elif decoy_strategy == "shifted reversal":
                    decoy_sequence = shift_reverse_sequence(sequence, shifted_amino_acids)

                decoy_locus_to_sequence_map[decoy_locus] = {'description':decoy_description, 'sequence':decoy_sequence}

            target_protein_col.metric(label="Target Proteins", value=len(locus_to_sequence_map))
            decoy_protein_col.metric(label="Decoy Proteins", value=len(decoy_locus_to_sequence_map))

            new_fasta_lines = fasta_from_locus_to_sequence_map({**locus_to_sequence_map, **decoy_locus_to_sequence_map})
            st.download_button("Download FASTA", "".join(new_fasta_lines), file_name=f"{Path(fasta_file.name).stem}_filter.fasta")

            st.balloons()
