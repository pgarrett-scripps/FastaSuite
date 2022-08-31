import re
import streamlit as st
import requests
from datetime import date

st.title("Download FASTA")
with st.expander("Help"):
     st.markdown("""
     FASTA files will be downloaded directly from uniprot. 
     
     **Input:**
     
     **Organisms**: The name of the organisms you would like to include in the in the Fasta file. 
     Selecting multiple will concatenate all proteins into a single FASTA file.
     
     **Reviewed**: If selected then the FASTA file will only contain reviewed proteins.
     """)

organims_to_tax_id = {
     'Human':9606,
     'E. coli':83333,
     'C. elegans':6239
}

def fasta_url_download_link(tax_id, is_reviewed):
     return f'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28taxonomy_id%3A{tax_id}%29%20AND%20%28reviewed%3A{str(is_reviewed).lower()}%29'

options = st.multiselect(
     'Organisms',
     list(organims_to_tax_id.keys())
)

is_reviewed = st.checkbox("reviewed",value=True)

if st.button("Get FASTA"):

     with st.spinner("Downloading fasta..."):

          responses = []
          for option in options:
               url = fasta_url_download_link(organims_to_tax_id[option], is_reviewed)
               response = requests.get(url)
               responses.append(response.content)

          combied_response = b'\n'.join(responses)

          today = date.today()
          mdy = today.strftime("%m_%d_%Y")
          fasta_file_name = '_'.join([re.sub(r'[^a-zA-Z]', '', option).lower() for option in options])
          fasta_file_name = f'{fasta_file_name}{"_reviewed" if is_reviewed else ""}_{mdy}.fasta'
          st.download_button(
               label=f"Download {fasta_file_name}",
               data=combied_response,
               file_name=fasta_file_name,
               mime='text/csv',
          )

