import copy
from collections import Counter
from pathlib import Path

from Bio.Seq import Seq
import streamlit as st
from Bio import SeqIO

from constants import VALID_AMINO_ACIDS, DECOY_GENERATOR_HELP_MESSAGE, DECOY_FLAG_HELP_MESSAGE, \
    RANDOM_SEED_HELP_MESSAGE, STATIC_AMINO_ACID_HELP_MESSAGE, MARKOV_STATE_SIZE_HELP_MESSAGE, KMER_SIZE_HELP_MESSAGE, \
    SHIFTED_AMINO_ACID_HELP_MESSAGE
from decoy_util import reverse_sequence, static_shuffle_sequence, shuffle_sequence, \
    exchange_sequence, shift_reverse_sequence, make_sequence_markov_model, make_locus_name_markov_model, \
    make_gene_name_markov_model, predict_sequence_from_markov_model, predict_locus_name_from_markov_model, \
    predict_gene_name_from_markov_model, construct_bruijn_graph, randomize_nodes, construct_sequence
from utils import map_locus_to_sequence_from_fasta, fasta_from_locus_to_sequence_map, write_temp_file, make_temp_file

st.title("Decoy Generator")
with st.expander("Help"):
    st.markdown(DECOY_GENERATOR_HELP_MESSAGE)

fasta_file = st.file_uploader(label="Upload FASTA", type=".fasta")

decoy_flag = st.text_input(label="Decoy Flag (set to 'Reverse_' for IP2)", value="DECOY_",
                           help=DECOY_FLAG_HELP_MESSAGE)
decoy_strategy = st.radio(label="Decoy Strategy",
                          options=('reverse', 'shuffle', 'markov', 'exchange', 'shifted reversal', 'deBruijn'),
                          index=0)

random_seed = None
if decoy_strategy in {'shuffle', 'markov', 'deBruijn'}:
    random_seed = st.number_input("Seed", value=123, help=RANDOM_SEED_HELP_MESSAGE)

static_amino_acids = None
if decoy_strategy in {'shuffle', 'deBruijn'}:
    static_amino_acids = st.multiselect(label='Static Residues', options=list(VALID_AMINO_ACIDS), default=None,
                                        help=STATIC_AMINO_ACID_HELP_MESSAGE)
markov_state_size = None
if decoy_strategy == 'markov':
    markov_state_size = st.number_input("Markov Chain Memory", min_value=2, max_value=2, value=2,
                                        help=MARKOV_STATE_SIZE_HELP_MESSAGE)

kmer_size = None
if decoy_strategy == 'deBruijn':
    kmer_size = st.number_input("K-mer size", min_value=2, max_value=5, value=2,
                                help=KMER_SIZE_HELP_MESSAGE)
aa_exchange_map = {}
if decoy_strategy == 'exchange':
    for aa in VALID_AMINO_ACIDS:
        aa_exchange_map[aa] = st.selectbox(label=f'Exchange {aa} with:', options=list(VALID_AMINO_ACIDS),
                                           index=list(VALID_AMINO_ACIDS).index(aa))

shifted_amino_acids = None
if decoy_strategy == 'shifted reversal':
    shifted_amino_acids = set(st.multiselect("Shifted Residues", options=list(VALID_AMINO_ACIDS),
                                             default=['K', 'R'], help=SHIFTED_AMINO_ACID_HELP_MESSAGE))

if st.button("Generate Decoys"):

    if not fasta_file:
        st.warning('Upload a FASTA file!')
        st.stop()

    file_type = None
    if fasta_file.name.endswith('.fasta'):
        file_type = 'fasta'

    elif fasta_file.name.endswith('.xml'):
        file_type = 'uniprot-xml'

    with write_temp_file(fasta_file) as tmp_file:
        with open(tmp_file) as handle:
            records = list(SeqIO.parse(handle, file_type))

    with st.spinner("Generating decoys..."):

        fasta_lines = fasta_file.getvalue().decode("utf-8").split("\n")

        sequences = [str(record.seq) for record in records]
        lengths = [len(sequence) for sequence in sequences]
        min_len, max_len = min(lengths), max(lengths)

        if decoy_strategy == 'markov':
            sequence_model = make_sequence_markov_model(sequences, markov_state_size)

        if decoy_strategy == 'deBruijn':
            nodes = construct_bruijn_graph(sequences, kmer_size)
            concatenated_sequence = "".join(sequences)
            amino_acids_count = Counter(concatenated_sequence)
            amino_acids_frequency = {aa: amino_acids_count[aa] / len(concatenated_sequence) for aa in
                                     amino_acids_count}
            randomize_nodes(nodes, static_amino_acids, amino_acids_frequency)

        decoy_records = []
        for record in records:

            sequence = str(record.seq)
            decoy_record = copy.deepcopy(record)
            decoy_record.name = decoy_flag + record.name
            decoy_record.id = decoy_flag + record.id

            if decoy_strategy == 'reverse':
                decoy_record.seq = Seq(reverse_sequence(sequence))
            elif decoy_strategy == 'shuffle':
                if static_amino_acids:
                    decoy_record.seq = Seq(static_shuffle_sequence(sequence, static_amino_acids))
                else:
                    decoy_record.seq = Seq(shuffle_sequence(sequence))
            elif decoy_strategy == 'markov':

                decoy_record.seq = Seq(predict_sequence_from_markov_model(sequence_model, min_len, max_len))

            elif decoy_strategy == 'exchange':
                decoy_record.seq = Seq(exchange_sequence(sequence, aa_exchange_map))

            elif decoy_strategy == "shifted reversal":
                decoy_record.seq = Seq(shift_reverse_sequence(sequence, shifted_amino_acids))

            elif decoy_strategy == 'deBruijn':
                decoy_record.seq = Seq(construct_sequence(nodes, sequence, kmer_size))

            decoy_records.append(decoy_record)

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

        with make_temp_file() as tmp_file:
            SeqIO.write(records + decoy_records, tmp_file, file_type)
            with open(tmp_file) as file:
                combined_contents = file.read()

        with make_temp_file() as tmp_file:
            SeqIO.write(decoy_records, tmp_file, file_type)
            with open(tmp_file) as file:
                decoy_contents = file.read()


        def create_fasta_file_name(combined: bool):
            return f"{Path(fasta_file.name).stem}_{decoy_strategy.replace(' ', '_')}" \
                          f"{param_tag}_combined_{combined}.fasta".lower()

        st.download_button(f"Download Combined",
                           combined_contents,
                           file_name=create_fasta_file_name(combined=True))

        st.download_button(f"Download Decoy",
                           decoy_contents,
                           file_name=create_fasta_file_name(combined=False))
