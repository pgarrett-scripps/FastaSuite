import streamlit as st
from pathlib import Path
from typing import Dict, List, Tuple
from collections import OrderedDict

def map_locus_to_sequence_from_fasta(fasta_lines):
    """Parse FASTA lines into locus-to-sequence mapping."""
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
    """Create FASTA string from locus-to-sequence mapping."""
    lines = []
    for locus in locus_to_sequence_map:
        lines.append(f'>{locus} {locus_to_sequence_map[locus]["description"]}\n')
        lines.append(f'{locus_to_sequence_map[locus]["sequence"]}\n')
    return ''.join(lines)

def get_lines_from_uploaded_file(file):
    """Extract lines from uploaded file."""
    return file.getvalue().decode("utf-8").split("\n")

def merge_fasta_files(fasta_files: List, 
                     dedup_method: str = "locus",
                     remove_contaminants: bool = False,
                     contaminant_prefix: str = "",
                     remove_reverse: bool = False,
                     reverse_prefix: str = "",
                     match_strings: List[str] = None) -> Tuple[Dict[str, Dict], Dict[str, int]]:
    """Merge multiple FASTA files with filtering."""
    
    all_sequences = OrderedDict()
    stats = {'total_input': 0, 'duplicates_removed': 0, 'contaminants_removed': 0, 
             'reverse_removed': 0, 'no_match_removed': 0, 'final_count': 0}
    
    seen_loci = set()
    seen_sequences = set()
    
    for fasta_file in fasta_files:
        fasta_lines = get_lines_from_uploaded_file(fasta_file)
        locus_to_sequence_map = map_locus_to_sequence_from_fasta(fasta_lines)
        
        for locus, sequence_data in locus_to_sequence_map.items():
            stats['total_input'] += 1
            
            # Filter by match strings first (if provided)
            if match_strings:
                has_match = any(match_str in locus for match_str in match_strings)
                if not has_match:
                    stats['no_match_removed'] += 1
                    continue
            
            # Filter contaminants
            if remove_contaminants and contaminant_prefix and contaminant_prefix in locus:
                stats['contaminants_removed'] += 1
                continue
                
            # Filter reverse sequences
            if remove_reverse and reverse_prefix and reverse_prefix in locus:
                stats['reverse_removed'] += 1
                continue
            
            # Check for duplicates
            is_duplicate = False
            if dedup_method == "locus":
                if locus in seen_loci:
                    is_duplicate = True
                else:
                    seen_loci.add(locus)
            elif dedup_method == "sequence":
                sequence = sequence_data['sequence']
                if sequence in seen_sequences:
                    is_duplicate = True
                else:
                    seen_sequences.add(sequence)
            
            if is_duplicate:
                stats['duplicates_removed'] += 1
                continue
            
            all_sequences[locus] = sequence_data
    
    stats['final_count'] = len(all_sequences)
    return all_sequences, stats

def main():
    st.title("FASTA Merger")
    
    # File uploads
    fasta_files = st.file_uploader("FASTA files", type=[".fasta", ".fa", ".fas", ".faa"], accept_multiple_files=True)
    
    # Match strings input
    st.subheader("Filter by protein IDs (optional)")
    match_input = st.text_area("Protein IDs to keep (one per line or space-separated)", 
                              placeholder="Q86XF7\nQ9UKM9\nP26022\nQ8WWC4")
    
    # Parse match strings
    match_strings = []
    if match_input.strip():
        # Split by newlines and spaces, filter empty strings
        raw_strings = match_input.replace('\n', ' ').split()
        match_strings = [s.strip() for s in raw_strings if s.strip()]
    
    if match_strings:
        st.info(f"Will keep proteins containing: {', '.join(match_strings[:10])}{'...' if len(match_strings) > 10 else ''}")
    
    # Options
    col1, col2 = st.columns(2)
    
    with col1:
        dedup_method = st.radio("Remove duplicates by:", ["locus", "sequence"])
        
    with col2:
        remove_contaminants = st.checkbox("Remove contaminants")
        contaminant_prefix = st.text_input("Contaminant prefix", "contaminant_", disabled=not remove_contaminants)
        
        remove_reverse = st.checkbox("Remove reverse/decoy")
        reverse_prefix = st.text_input("Reverse prefix", "Reverse_", disabled=not remove_reverse)
    
    # Process
    if st.button("Merge"):
        if not fasta_files:
            st.error("Upload FASTA files")
            return
        
        merged_sequences, stats = merge_fasta_files(
            fasta_files=fasta_files,
            dedup_method=dedup_method,
            remove_contaminants=remove_contaminants,
            contaminant_prefix=contaminant_prefix if remove_contaminants else "",
            remove_reverse=remove_reverse,
            reverse_prefix=reverse_prefix if remove_reverse else "",
            match_strings=match_strings
        )
        
        # Results
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Input", stats['total_input'])
        col2.metric("Duplicates", stats['duplicates_removed'])
        col3.metric("Filtered", stats['contaminants_removed'] + stats['reverse_removed'] + stats['no_match_removed'])
        col4.metric("Final", stats['final_count'])
        
        if merged_sequences:
            merged_fasta = fasta_from_locus_to_sequence_map(merged_sequences)
            st.download_button("Download", merged_fasta, "merged.fasta")
        else:
            st.warning("No sequences remaining")

if __name__ == "__main__":
    main()