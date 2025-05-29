from Bio import SeqIO
import glob
import re


def merge_and_update_gff(file_list, output_file):
    with open(output_file, 'w', encoding='utf-8') as outfile:
        for file in sorted(file_list):
            with open(file, 'r', encoding='utf-8') as infile:
                outfile.write(f"# Source: {file}\n")
                for line in infile:
                    if line.strip() and not line.startswith('#'):
                        seq_id = line.split('\t')[0]
                        match = re.match(r'(SEQF\d+)', seq_id)
                        prefix = match.group(1) if match else None
                        if prefix:
                            updated_line = re.sub(
                                r'(locus_tag=)(pgaptmp_\d+)',
                                rf'\1{prefix}-\2',
                                line
                            )
                            # update ID
                            updated_line = re.sub(
                                r'(ID=)((gene|cds)-pgaptmp_\d+)',
                                rf'\1{prefix}-\2',
                                updated_line
                            )
                            # update Parent
                            updated_line = re.sub(
                                r'(Parent=)(gene-pgaptmp_\d+)',
                                rf'\1{prefix}-\2',
                                updated_line
                            )
                            outfile.write(updated_line)
                        else:
                            outfile.write(line)
gff_files = glob.glob('./*_modified_annot.gff')

output_file = "combined.gff"

merge_and_update_gff(gff_files, output_file)
