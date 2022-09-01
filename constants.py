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

# FASTA SCORE
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

# FASTA DECOY
DECOY_GENERATOR_HELP_MESSAGE = """
    Generate Decoy Proteins in a FASTA file

    **Input**

    **FASTA File**: FASTA file containing proteins. All proteins in this file will be used to generate a decoy 
    variant so there should be no decoy sequences present
    
    **Decoy Flag:** The decoy flag to prepend locus names ("Reverse_" -> ">Reverse_sp|XXXX|YYYY")
    Set to "Reverse_" for IP2
    
    **Decoy Strategies**
    
    - **reverse:** reverse the protein sequence
    - **shuffle:** randomly shuffle the protein sequence
    - **shuffle static:** randomly shuffle the residues between given static residues
    - **markov:** trains a markov chain model on target proteins and uses this model to predict decoy proteins (Think autocompletion)
    - **exchange:** exchange residues with their given replacement
    - **shifted reversal:** reverse the sequence, then switch the given amino acids with their predecessor 
    - **DeBruijn:** randomly shuffle amino acids while keeping repeated patterns (similar to shuffle but repeat sequences are preserved)
    - **DeBruijn static:** randomly shuffle amino acids while keeping repeated patterns and static residues (similar to static shuffle but repeat sequences are preserved)
    
    
    The "ideal" decoy database should statistically mimic the target database. Such a database should conserve 
    amino acid frequency, peptide lengths, peptide masses, protein lengths, repeat sequences within proteins, and 
    repeat sequences between proteins.
    """

RANDOM_SEED_HELP_MESSAGE = "Used to initialize a pseudorandom number generator. Keeping the same number will allow for reproducible decoy generation."
STATIC_AMINO_ACID_HELP_MESSAGE = "Static Residues will have the same sequential location in decoy and target proteins"
MARKOV_STATE_SIZE_HELP_MESSAGE = 'The memory of the markov chain: how many previous residues to consider when predicting next amino acid'
KMER_SIZE_HELP_MESSAGE= 'The number of residues needed to identify repeat sequences. The first N elements in a repeated sequence will lost.'
SHIFTED_AMINO_ACID_HELP_MESSAGE = "These residues will be swapped with the next residue. Then the sequence will be reversed."