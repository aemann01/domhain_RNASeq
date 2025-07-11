from Bio import SeqIO
import glob
import re
import os
from io import StringIO

def clean_locus_lines(text):
    cleaned_lines = []
    for line in text.splitlines(keepends=True):
        if line.startswith("LOCUS"):
            # Remove any unexpected characters (e.g. '>')
            line = re.sub(r'[^\w\s\-]', '', line)
        cleaned_lines.append(line)
    return ''.join(cleaned_lines)

def add_molecule_type(record):
    # Check if 'molecule_type' is missing, and add a default if so
    if 'molecule_type' not in record.annotations:
        record.annotations['molecule_type'] = 'DNA' 
    return record

def combine_genbank_files(input_files, output_file):
    records = []

    for input_file in input_files:
        print(f"Processing {input_file}...")
        try:
            with open(input_file, 'r') as infile:
                raw_text = infile.read()
                cleaned_text = clean_locus_lines(raw_text)
                handle = StringIO(cleaned_text)
                parsed = list(SeqIO.parse(handle, 'genbank'))
                
                if parsed:
                    # Add molecule_type annotation to each record
                    for record in parsed:
                        record = add_molecule_type(record)
                    records.extend(parsed)
                else:
                    print(f"No valid records in: {input_file}")
        except Exception as e:
            print(f"❌ Failed to parse {input_file}: {e}")

    if records:
        with open(output_file, 'w') as outfile:
            SeqIO.write(records, outfile, 'genbank')
        print(f"Combined GenBank file saved as: {output_file}")
    else:
        print("No records to write. Check your input files.")

# Get list of GenBank files
gbk_files = glob.glob(os.path.expanduser('~/rna_dohmain/11-perio/09-GO_terms/annotations/*_modified_annot.gbk'))

output_file = 'combined_homd.gbk'
combine_genbank_files(gbk_files, output_file)
