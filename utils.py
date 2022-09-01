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
    return f'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28taxonomy_id%3A{tax_id}%29%20AND%20%28reviewed%3A{str(is_reviewed).lower()}%29'


