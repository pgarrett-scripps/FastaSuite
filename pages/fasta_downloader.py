import streamlit as st
import requests

st.subheader("Download FASTA file")
st.write("FASTA files will be downloaded directly from uniprot")

organims_to_tax_id = {
     'human':9606,
     'ecoli':83333,
     'celegans':6239
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

          st.download_button(
               label="Download Fasta",
               data=combied_response,
               file_name=f'{"_".join(options)}_reviewed_{str(is_reviewed).lower()}.fasta',
               mime='text/csv',
          )

