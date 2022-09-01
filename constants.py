VALID_AMINO_ACIDS = 'ARNDCEQGHILKMFPSTWYVUXZBJO'

# FASTA FILTER
FASTA_FILTERED_HELP_MESSAGE = """
    This app will generate a FASTA file containing only the proteins identified in the DTASelect-filter.txt file. 
    You can then upload the filtered FASTA file to IP2 to use a Library Search.
    
    **Input**
    
    **FASTA File**: Fasta file to filter (should be identical to the one used in the experiment)
    
    **Filter Files**: The files which contain identified proteins to filter the FASTA file (you can upload multiple)
    
    **Decoy & Contaminant Flags**: The flag which identifies decoy/contaminant proteins in the uploaded files. Protiens which contain this flag will be excluded from the filtered FASTA. ("Reverse_" for ">Reverse_sp|XXX|YYY" or "contaminant_" for ">contaminant_sp|XXX|YYY")
    
    **Note**: The DTASelect-filter.txt can contain Decoy and Contaminant proteins. 
    
    - Remove Decoy/Contaminant proteins from the filtered FASTA: set the flag inputs to the correct prefix.
    
    - Include Decoy/Contaminant proteins in filtered FASTA: leave the flag inputs blank.
    """

DECOY_FLAG_HELP_MESSAGE = "Flag used to identify Decoy peptides (if fasta contains '>DECOY_sp|XXXX|YYYY' set to 'DECOY_')"
CONTAMINANT_FLAG_HELP_MESSAGE = "Flag used to identify Contaminant peptides (if fasta contains '>contaminant_sp|XXXX|YYYY' set to 'contaminant_')"

# FASTA DOWNLOADER
FASTA_DOWNLOADER_HELP_MESSAGE= """
     FASTA files will be downloaded directly from uniprot. 

     **Input:**

     **Organisms**: The name of the organisms you would like to include in the in the Fasta file. 
     Selecting multiple will concatenate all proteins into a single FASTA file.

     **Reviewed**: If selected then the FASTA file will only contain reviewed proteins.
     """

ORGANISM_TO_TAXONOMY_ID = {
     'Human': 9606,
     'E. coli': 83333,
     'C. elegans': 6239
}

FASTA_SCORE_HELP_MESSAGE = """
        Scores a FASTA file according to decoy/target statistics.
        The "ideal" FASTA file contains target and decoy proteins which statistically resemble each other.

        **Input:**

        **Enzyme Sites**: The digestion sites. Use K,R for trypsin. 

        **Min/Max Peptide Lengths**: Keep only digested peptides within these bounds (bounds are inclusive)

        **Decoy Flag**: The flag used to represent Decoy Proteins
    """
ENZYME_HELP_MESSAGE = "Residues to cleave after. For Trypsin use (K & R), for nonspecific leave blank"
MIN_MAX_PEPTIDE_LENGTH_HELP_MESSAGE = "Sets range of supported peptides (min <= len_of_peptide <= max)"
