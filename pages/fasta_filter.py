from pathlib import Path

import streamlit as st
from senpy.dtaSelectFilter.parser2 import read_file

from utils import map_locus_to_sequence_from_fasta, fasta_from_locus_to_sequence_map

st.title("Filter FASTA")
with st.expander("Help"):
    st.markdown("""
    
    This app will generate a FASTA file containing only identified proteins. You can then upload the filtered FASTA file  to IP2 to use a Library Search.
    
    **Input**
    
    **FASTA File**: Fasta file to filter (should be identical to the one used in the experiment)
    
    **Filter Files**: The files which contain identified proteins to filter the FASTA file (you can upload multiple)
    
    **Decoy & Contaminant Flags**: The flag which identifies decoy/contaminant proteins in the uploaded files. Protiens which contain this flag will be excluded from the filtered FASTA. ("Reverse_" for ">Reverse_sp|XXX|YYY" or "contaminant_" for ">contaminant_sp|XXX|YYY")
    
    *To Include Decoy or Contaminant proteins: leave the Flag input blank.
    """)

fasta_file = st.file_uploader("FASTA", ".fasta")
dta_filter_files = st.file_uploader("DTASelect-filter.txt", ".txt", accept_multiple_files=True)

decoy_flag = st.text_input("Decoy Flag", "Reverse_", help="Flag used to identify Decoy peptides (if fasta contains"
                                                          " '>DECOY_sp|XXXX|YYYY' set to 'DECOY_')")
contaminant_flag = st.text_input("Contaminant Flag", "contaminant_",
                                 help="Flag used to identify Contaminant peptides (if fasta contains "
                                      "'>contaminant_sp|XXXX|YYYY' set to 'contaminant_')")

if st.button("Generate"):

    if not fasta_file:
        st.warning("Upload FASTA file...")
    if not dta_filter_files:
        st.warning("Upload DTASelect-filter.txt file...")

    if fasta_file and dta_filter_files:

        fasta_lines = fasta_file.getvalue().decode("utf-8").split("\n")
        results_list = []
        for dta_filter_file in dta_filter_files:
            filter_lines = dta_filter_file.getvalue().decode("utf-8").split("\n")
            _, results, _ = read_file(filter_lines)
            results_list.extend(results)

        locus_to_sequence_map = map_locus_to_sequence_from_fasta(fasta_lines)
        protein_locuses = list({protein.locus for result in results_list for protein in result.proteins})
        if decoy_flag and decoy_flag != "":
            protein_locuses = [locus for locus in protein_locuses if decoy_flag not in locus]
        if contaminant_flag and contaminant_flag != "":
            protein_locuses = [locus for locus in protein_locuses if contaminant_flag not in locus]

        dta_filter_locus_to_sequence_map = {locus: locus_to_sequence_map[locus] for locus in protein_locuses}

        new_fasta_lines = fasta_from_locus_to_sequence_map(dta_filter_locus_to_sequence_map)
        st.download_button("Download FASTA",
                           "".join(new_fasta_lines),
                           file_name=f"{Path(fasta_file.name).stem}_filter.fasta")
