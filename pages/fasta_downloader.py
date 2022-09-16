import re
import streamlit as st
import requests
from datetime import date

from constants import FASTA_DOWNLOADER_HELP_MESSAGE, ORGANISM_TO_TAXONOMY_ID
from utils import fasta_url_download_link

st.title("Download FASTA")
with st.expander("Help"):
    st.markdown(FASTA_DOWNLOADER_HELP_MESSAGE)

options = st.multiselect(
    label='Organisms',
    options=list(ORGANISM_TO_TAXONOMY_ID)
)

is_reviewed = st.checkbox("reviewed", value=True)

if st.button("Get FASTA"):

    with st.spinner("Downloading fasta..."):

        responses = []
        for option in options:
            url = fasta_url_download_link(ORGANISM_TO_TAXONOMY_ID[option], is_reviewed)
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
