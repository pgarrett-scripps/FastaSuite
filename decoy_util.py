import random
import re

import markovify

VALID_AMINO_ACIDS = 'ARNDCEQGHILKMFPSTWYVUXZBJO'


def reverse_sequence(sequence):
    return sequence[::-1]


def shuffle_sequence(sequence):
    return ''.join(random.sample(sequence, len(sequence)))


def static_shuffle_sequence(sequence, static_amino_acids=None):
    if static_amino_acids is None or len(static_amino_acids) == 0:
        return shuffle_sequence(sequence)

    aa_to_indexes = {aa: [m.start() for m in re.finditer(aa, sequence)] for aa in static_amino_acids}
    aa_indexes = [(index, aa) for aa in aa_to_indexes for index in aa_to_indexes[aa]]
    aas = [aa for i, aa in sorted(aa_indexes, key=lambda pair: pair[0])]
    sub_sequences = re.split("|".join(static_amino_acids), sequence)
    decoy_sub_sequences = ["".join(random.sample(sub_sequence, len(sub_sequence))) for sub_sequence in sub_sequences]

    for i in reversed(range(1, len(decoy_sub_sequences))):
        decoy_sub_sequences.insert(i, aas[i - 1])

    return ''.join(decoy_sub_sequences)


def exchange_sequence(sequence, amino_acid_exchange_map=None):
    if amino_acid_exchange_map is None or amino_acid_exchange_map == {}:
        return sequence

    return "".join([amino_acid_exchange_map.get(aa, aa) for aa in sequence])


def shift_reverse_sequence(sequence, shift_amino_acids=None):

    if shift_amino_acids is None or len(shift_amino_acids) == 0:
       return sequence

    decoy_sequence = list(reverse_sequence(sequence))
    prev_aa = decoy_sequence[-1]
    for i in range(len(decoy_sequence)):
        aa = decoy_sequence[i]
        if aa in shift_amino_acids:
            if i > 0:
                decoy_sequence[i-1] = decoy_sequence[i]
                decoy_sequence[i] = prev_aa
        prev_aa = aa
    return "".join(decoy_sequence)


def make_sequence_markov_model(sequences, markov_state_size):
    protein_sentences = [" ".join(list(sequence)) + "\n" for sequence in sequences]
    return markovify.NewlineText("".join(protein_sentences), state_size=markov_state_size)


def make_locus_name_markov_model(locus_names):
    locus_name_sentences = [" ".join(list(locus_name))+"\n" for locus_name in locus_names]
    return markovify.NewlineText("".join(locus_name_sentences))


def make_gene_name_markov_model(gene_names):
    gene_name_sentences = [" ".join(list(gene_name)) + "\n" for gene_name in gene_names]
    return markovify.NewlineText("".join(gene_name_sentences))


def predict_sequence_from_markov_model(model, min_len, max_len):
    return model.make_sentence(min_chars=min_len, max_chars=max_len, test_output=False).replace(" ", "")


def predict_locus_name_from_markov_model(model, min_chars, max_chars):
    return model.make_sentence(min_chars=min_chars, max_chars=max_chars, test_output=False).replace(" ", "")


def predict_gene_name_from_markov_model(model, min_chars, max_chars):
    return model.make_sentence(min_chars=min_chars, max_chars=max_chars, test_output=False).replace(" ", "")


def build_kmers(sequence, k):
    kmers = []
    n_kmers = len(sequence) - k + 1
    for i in range(n_kmers):
        kmer = sequence[i:i + k]
        kmers.append(kmer)
    return kmers


def construct_bruijn_graph(sequences, k):
    nodes = {}

    for sequence in sequences:
        sequence = "".join(['-'] * k) + sequence

        kmers = build_kmers(sequence, k)
        for i, kmer in enumerate(kmers):
            if i == len(kmers) - 1:
                break
            next_kmer = kmers[i + 1]
            next_letter = next_kmer[-1]
            if kmer not in nodes:
                nodes[kmer] = {next_kmer: next_letter}
            else:
                nodes[kmer][next_kmer] = next_letter
    return nodes


def construct_sequence(nodes, sequence, k):
    new_sequence = []

    sequence = "".join(['-'] * k) + sequence
    kmers = build_kmers(sequence, k)
    for i, kmer in enumerate(kmers):
        if i == len(kmers) - 1:
            break
        next_kmer = kmers[i + 1]
        next_letter = next_kmer[-1]

        letter = nodes[kmer][next_kmer]
        new_sequence.append(letter)

    return "".join(new_sequence)


def randomize_nodes(nodes, static_residues=None, amino_acids_frequency=None):
    amino_acid_weights = [amino_acids_frequency.get(aa, 0) for aa in VALID_AMINO_ACIDS]
    for aa in static_residues:
        amino_acid_weights[VALID_AMINO_ACIDS.index(aa)] = 0

    for key1 in nodes:
        for key2 in nodes[key1]:
            aa = nodes[key1][key2]
            if aa in static_residues:
                continue
            nodes[key1][key2] = random.choices(VALID_AMINO_ACIDS, weights=amino_acid_weights, k=1)[0]






