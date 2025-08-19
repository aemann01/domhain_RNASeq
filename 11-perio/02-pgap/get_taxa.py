import csv

def get_seqid_info(annots_file, seqid_info_file, output_file):
    # make it a dictionary to speed up 
    seqid_dict = {}

    # read the seqid_info_file only once and build the seqid_dict
    with open(seqid_info_file, mode='r') as seqid_file:
        csv_reader = csv.DictReader(seqid_file)
        for row in csv_reader:
            seq_id = row['SEQ_ID'].split('_')[0] 
            seqid_dict[seq_id] = row

    # open the annots file and output file
    with open(annots_file, mode='r') as annots_file, open(output_file, mode='w', newline='') as output:
        annots_lines = annots_file.readlines()
        writer = csv.writer(output)

        header = ['SEQ_ID', 'HMT_ID', 'Genus', 'Species', 'Strain', 'Contigs', 'Combined_Size', 'Sequence_Source']
        writer.writerow(header)

        # process each line in annots.txt
        output_lines = []
        for line in annots_lines:
            seqid_part = line.split('_')[0]  # part before first underscore

            # Check for exact match in seqid_dict first
            if seqid_part in seqid_dict:
                output_lines.append(list(seqid_dict[seqid_part].values()))  # append matched row
            else:
                # If no exact match, check for a partial match
                match_found = False
                for row in seqid_dict.values():
                    if seqid_part in row['SEQ_ID']:  # check for partial match in SEQ_ID
                        output_lines.append(list(row.values()))  # append matched row
                        match_found = True
                        break

                # if no match
                if not match_found:
                    output_lines.append([seqid_part, 'Not Found'] + ['Not Found'] * (len(header) - 1))

        writer.writerows(output_lines)

if __name__ == '__main__':
    annots_file = 'annots.txt'
    seqid_info_file = 'SEQID_info.csv'
    output_file = 'taxa.csv'

    get_seqid_info(annots_file, seqid_info_file, output_file)

    print(f"Output written to {output_file}")
