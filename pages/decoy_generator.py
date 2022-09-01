from collections import Counter
from pathlib import Path

import streamlit as st

from constants import VALID_AMINO_ACIDS
from decoy_util import reverse_sequence, static_shuffle_sequence, shuffle_sequence, \
    exchange_sequence, shift_reverse_sequence, make_sequence_markov_model, make_locus_name_markov_model, \
    make_gene_name_markov_model, predict_sequence_from_markov_model, predict_locus_name_from_markov_model, \
    predict_gene_name_from_markov_model, construct_bruijn_graph, randomize_nodes, construct_sequence
from utils import map_locus_to_sequence_from_fasta, fasta_from_locus_to_sequence_map

st.title("Decoy Generator")
with st.expander("Help"):
    st.markdown("""

    Generate Decoy Proteins in a FASTA file

    **Input**

    **FASTA File**: FASTA file containing proteins. All proteins in this file will be used to generate a decoy 
    variant so there should be no decoy sequences present
    
    **Decoy Flag:** The decoy flag to prepend locus names ("Reverse_" -> ">Reverse_sp|XXXX|YYYY")
    Set to "Reverse_" for IP2
    
    **Decoy Strategies**
    
    - **reverse:** reverse the sequence
    - **shuffle:** randomly shuffle the sequence
    - **shuffle static:** randomly shuffle residues between given static residues
    - **markov:** trains a markov chain model on target proteins and uses this model to predict decoy proteins
    - **exchange:** exchange the provided residues with the replacement
    - **shifted reversal:** reverse the sequence, then switch the given amino acids with their predecessor 
    - **DeBruijn:** randomly shuffle amino acids while keeping repeated patterns 
    - **DeBruijn static:** randomly shuffle amino acids while keeping repeated patterns and static residues


    
    The "ideal" decoy database should statistically mimic the target database. Such a database should conserve 
    amino acid frequency, peptide lengths, peptide masses, protein lengths, repeat sequences within proteins, and 
    repeat sequences between proteins.

    """)

fasta_file = st.file_uploader("Upload FASTA", type=".fasta")

decoy_flag = st.text_input("Decoy Flag (set to 'Reverse_' for IP2)", "DECOY_",
                           help="The flag used to identify Decoy peptides: 'DECOY_' "
                                "will result in '>DECOY_sp|XXXX|YYYY'")
decoy_strategy = st.radio("Decoy Strategy",
                          ('reverse', 'shuffle', 'markov', 'exchange', 'shifted reversal', 'deBruijn'), index=0)

random_seed = None
if decoy_strategy in {'shuffle', 'markov', 'deBruijn'}:
    random_seed = st.number_input("Seed", value=123,
                                  help="Used to initialize a pseudorandom number generator. Keeping the same "
                                       "number will allow for reproducible decoy generation.")

static_amino_acids = None
if decoy_strategy in {'shuffle', 'deBruijn'}:
    static_amino_acids = st.multiselect(label='Static Residues', options=list(VALID_AMINO_ACIDS), default=None,
                                        help="Static Residues will have the same sequential location in decoy and"
                                             " target proteins")
markov_state_size = None
if decoy_strategy == 'markov':
    markov_state_size = st.number_input("Markov Chain Memory", min_value=2, max_value=2, value=2,
                                        help='The memory of the markov chain: how many previous residues to consider '
                                             'when predicting next amino acid')

kmer_size = None
if decoy_strategy == 'deBruijn':
    kmer_size = st.number_input("K-mer size", min_value=2, max_value=5, value=2,
                                help='The number of residues needed to identify repeat sequences. '
                                     'The first N elements in a repeated sequence will lost.')
aa_exchange_map = {}
if decoy_strategy == 'exchange':
    for aa in VALID_AMINO_ACIDS:
        aa_exchange_map[aa] = st.selectbox(label=f'Exchange {aa} with:', options=list(VALID_AMINO_ACIDS),
                                           index=list(VALID_AMINO_ACIDS).index(aa))

shifted_amino_acids = None
if decoy_strategy == 'shifted reversal':
    shifted_amino_acids = set(st.multiselect("Shifted Residues", options=list(VALID_AMINO_ACIDS),
                                             default=['K', 'R'],
                                             help="These residues will be swapped with the next residue. "
                                                  "Then the sequence will be reversed."))

if st.button("Generate Decoys"):

    if not fasta_file:
        st.warning('Upload a FASTA file!')

    if fasta_file is not None:

        with st.spinner("Generating decoys..."):

            fasta_lines = fasta_file.getvalue().decode("utf-8").split("\n")

            locus_to_sequence_map = map_locus_to_sequence_from_fasta(fasta_lines)
            lengths = [len(locus_to_sequence_map[locus]['sequence']) for locus in locus_to_sequence_map]
            min_len, max_len = min(lengths), max(lengths)
            org_name = list(locus_to_sequence_map.keys())[0].split("|")[2].split("_")[1]
            sequences = [locus_to_sequence_map[locus]['sequence'] for locus in locus_to_sequence_map]

            if decoy_strategy == 'markov':
                sequence_model = make_sequence_markov_model(sequences, markov_state_size)

                locus_names = [locus.split("|")[1] for locus in locus_to_sequence_map]
                locus_name_model = make_locus_name_markov_model(locus_names)

                gene_names = [locus.split("|")[1].split("_")[0] for locus in locus_to_sequence_map]
                gene_name_model = make_gene_name_markov_model(gene_names)

            if decoy_strategy == 'deBruijn':
                nodes = construct_bruijn_graph(sequences, kmer_size)
                concatenated_sequence = "".join(sequences)
                amino_acids_count = Counter(concatenated_sequence)
                amino_acids_frequency = {aa: amino_acids_count[aa] / len(concatenated_sequence) for aa in
                                         amino_acids_count}
                randomize_nodes(nodes, static_amino_acids, amino_acids_frequency)

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

                elif decoy_strategy == 'deBruijn':
                    decoy_sequence = construct_sequence(nodes, sequence, kmer_size)

                decoy_locus_to_sequence_map[decoy_locus] = {'description': decoy_description,
                                                            'sequence': decoy_sequence}

            static_tag = ""
            if static_amino_acids:
                aa_tag = ''.join([aa for aa in static_amino_acids])
                static_tag = f"_static_{aa_tag}"

            shifted_tag = ""
            if shifted_amino_acids:
                aa_tag = ''.join([aa for aa in shifted_amino_acids])
                shifted_tag = f"_{aa_tag}"

            exchange_tag = ""
            if aa_exchange_map:
                aa_tag = '_'.join([f"{aa}{aa_exchange_map[aa]}" for aa in aa_exchange_map if aa != aa_exchange_map[aa]])
                exchange_tag = f"_{aa_tag}"

            markov_memory_tag = ""
            if markov_state_size:
                markov_memory_tag = f'_memory_{markov_state_size}'

            random_seed_tag = ""
            if random_seed:
                random_seed_tag = f'_seed_{random_seed}'

            kmer_size_tag = ""
            if kmer_size:
                kmer_size_tag = f'_kmer_{kmer_size}'

            param_tag = f'{static_tag}{shifted_tag}{exchange_tag}{markov_memory_tag}{random_seed_tag}{kmer_size_tag}'
            new_fasta_lines = fasta_from_locus_to_sequence_map({**locus_to_sequence_map, **decoy_locus_to_sequence_map})
            fasta_file_name = f"{Path(fasta_file.name).stem}_{decoy_strategy.replace(' ', '_')}{param_tag}.fasta".lower()
            st.download_button(f"Download {fasta_file_name}",
                               "".join(new_fasta_lines),
                               file_name=fasta_file_name)

            st.balloons()
