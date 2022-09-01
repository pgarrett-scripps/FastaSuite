from pathlib import Path

import streamlit as st
from senpy.dtaSelectFilter.parser2 import read_file

from constants import FASTA_FILTERED_HELP_MESSAGE, DECOY_FLAG_HELP_MESSAGE, CONTAMINANT_FLAG_HELP_MESSAGE
from utils import map_locus_to_sequence_from_fasta, fasta_from_locus_to_sequence_map, get_lines_from_uploaded_file

st.title("Filter FASTA")
with st.expander("Help"):
    st.markdown(FASTA_FILTERED_HELP_MESSAGE)

# File Input
fasta_file = st.file_uploader(label="FASTA", type=".fasta")
dta_filter_files = st.file_uploader(label="DTASelect-filter.txt", type=".txt", accept_multiple_files=True)

# Flag Input
decoy_flag = st.text_input(label="Decoy Flag", value="Reverse_", help=DECOY_FLAG_HELP_MESSAGE)
contaminant_flag = st.text_input(label="Contaminant Flag", value="contaminant_", help=CONTAMINANT_FLAG_HELP_MESSAGE)


if st.button("Filter"):

    if not fasta_file:
        st.warning("Upload FASTA file...")
    if not dta_filter_files:
        st.warning("Upload DTASelect-filter.txt file...")

    if fasta_file and dta_filter_files:

        # Parse fasta file
        fasta_lines = get_lines_from_uploaded_file(fasta_file)
        locus_to_sequence_map = map_locus_to_sequence_from_fasta(fasta_lines)

        # Parse DTASelect-filter and get a set of all protein locuses
        protein_locuses = set()
        for dta_filter_file in dta_filter_files:
            dta_filter_lines = get_lines_from_uploaded_file(dta_filter_file)
            _, dta_filter_results, _ = read_file(dta_filter_lines)
            for dta_filter_result in dta_filter_results:
                protein_locuses.update({protein_line.locus for protein_line in dta_filter_result.proteins})

        # Remove locuses matching these decoya nd contaminant flag.
        if decoy_flag: # *Note python interprets an empty string as falsy
            protein_locuses = [locus for locus in protein_locuses if decoy_flag not in locus]
        if contaminant_flag: # *Note python interprets an empty string as falsy
            protein_locuses = [locus for locus in protein_locuses if contaminant_flag not in locus]

        dta_filter_locus_to_sequence_map = {locus: locus_to_sequence_map[locus] for locus in protein_locuses}
        filtered_fasta = fasta_from_locus_to_sequence_map(dta_filter_locus_to_sequence_map)
        filtered_fasta_file_name = f"{Path(fasta_file.name).stem}_filter.fasta"
        st.download_button(label="Download FASTA", data=filtered_fasta, file_name=filtered_fasta_file_name)
