#!/usr/bin/env python3
import sys

def find_section(content, category):
    #Find a section in AGAT txt output
    lines = content.split('\n')
    section_lines = []
    in_section = False
    
    for line in lines:
        if '-' * 20 + f' {category} ' + '-' * 20 in line:
            in_section = True
            continue
        elif in_section and line.startswith('--------------------'):
            break
        elif in_section and line.strip():
            section_lines.append(line)
    
    return section_lines

def get_value(lines, key):
    #Get value from a line starting with the given key
    for line in lines:
        if line.strip().startswith(key):
            return line.strip().split()[-1]
    return '0'

def main(filename):
    # Read the input file
    with open(filename, 'r') as f:
        content = f.read()
    
    # Print header
    print("Category\tNo. genes\tNo. transcripts\tMean gene length (bp)\tNo. single-exon genes\tMean exons per transcript")
    
    # Process mrna (Protein-coding)
    mrna = find_section(content, 'mrna')
    if mrna:
        genes = get_value(mrna, 'Number of gene')
        transcripts = get_value(mrna, 'Number of mrna')
        length = get_value(mrna, 'mean gene length (bp)')
        single_exon = get_value(mrna, 'Number of single exon gene')
        mean_exons = get_value(mrna, 'mean exons per mrna')
        print(f"Protein-coding\t{genes}\t{transcripts}\t{length}\t{single_exon}\t{mean_exons}")
    
    # Process Ig/TCR segments (sum of c, d, j, v segments)
    segments = ['c_gene_segment', 'd_gene_segment', 'j_gene_segment', 'v_gene_segment']
    total_genes = 0
    total_transcripts = 0
    total_single_exon = 0
    lengths = []
    mean_exons = []
    
    for seg in segments:
        section = find_section(content, seg)
        if section:
            total_genes += int(get_value(section, 'Number of gene'))
            total_transcripts += int(get_value(section, f'Number of {seg}'))
            total_single_exon += int(get_value(section, 'Number of single exon gene'))
            length = float(get_value(section, 'mean gene length (bp)'))
            exons = float(get_value(section, f'mean exons per {seg}'))
            if length > 0:
                lengths.append(length)
            if exons > 0:
                mean_exons.append(exons)
    
    length_range = f"{int(min(lengths))}-{int(max(lengths))}" if lengths else "0"
    exon_range = f"{min(mean_exons):.1f}-{max(mean_exons):.1f}" if mean_exons else "0"
    print(f"Ig/TCR segments\t{total_genes}\t{total_transcripts}\t{length_range}\t{total_single_exon}\t{exon_range}")
    
    # Process pseudogenes
    pseudo = find_section(content, 'pseudogenic_transcript')
    if pseudo:
        genes = get_value(pseudo, 'Number of pseudogene')
        transcripts = get_value(pseudo, 'Number of pseudogenic_transcript')
        length = get_value(pseudo, 'mean pseudogene length (bp)')
        single_exon = get_value(pseudo, 'Number of single exon pseudogene')
        mean_exons = get_value(pseudo, 'mean exons per pseudogenic_transcript')
        print(f"Pseudogenes\t{genes}\t{transcripts}\t{length}\t{single_exon}\t{mean_exons}")
    
    # Process non-coding RNAs
    ncrna_types = {
        'lncRNA': 'lnc_rna',
        'snRNA': 'snrna',
        'snoRNA': 'snorna',
        'rRNA': 'rrna',
        'tRNA': 'trna',
        'miRNA': 'mirna',
        'scRNA': 'scrna'
    }
    
    for display_name, section_name in ncrna_types.items():
        section = find_section(content, section_name)
        if section:
            genes = get_value(section, 'Number of ncrna_gene')
            transcripts = get_value(section, f'Number of {section_name}')
            length = get_value(section, 'mean ncrna_gene length (bp)')
            single_exon = get_value(section, 'Number of single exon ncrna_gene')
            mean_exons = get_value(section, f'mean exons per {section_name}')
            print(f"{display_name}\t{genes}\t{transcripts}\t{length}\t{single_exon}\t{mean_exons}")
    
    # Process other non-coding (RNA + transcript)
    other_types = ['rna', 'transcript']
    total_genes = 0
    total_transcripts = 0
    total_single_exon = 0
    lengths = []
    mean_exons = []
    
    for type_name in other_types:
        section = find_section(content, type_name)
        if section:
            total_genes += int(get_value(section, 'Number of gene' if type_name == 'rna' else 'Number of ncrna_gene'))
            total_transcripts += int(get_value(section, f'Number of {type_name}'))
            total_single_exon += int(get_value(section, 'Number of single exon gene' if type_name == 'rna' else 'Number of single exon ncrna_gene'))
            length = float(get_value(section, 'mean gene length (bp)' if type_name == 'rna' else 'mean ncrna_gene length (bp)'))
            exons = float(get_value(section, f'mean exons per {type_name}'))
            if length > 0:
                lengths.append(length)
            if exons > 0:
                mean_exons.append(exons)
    
    length_range = f"{int(min(lengths))}-{int(max(lengths))}" if lengths else "0"
    exon_range = f"{min(mean_exons):.1f}-{max(mean_exons):.1f}" if mean_exons else "0"
    print(f"Other non-coding\t{total_genes}\t{total_transcripts}\t{length_range}\t{total_single_exon}\t{exon_range}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python make_agat_table.py <input_file>")
        sys.exit(1)
    main(sys.argv[1])
