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
    return lines