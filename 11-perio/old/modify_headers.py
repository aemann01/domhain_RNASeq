import os
import csv
# read
info_file = 'SEQID_info.csv'

# get the seqod into a dictionary
seq_info = {}
with open(info_file, 'r') as f:
    reader = csv.reader(f)
    next(reader)  
    for row in reader:
        seq_id = row[0].strip('"')
        genus = row[2].strip('"')
        species = row[3].strip('"')
        seq_info[seq_id] = (genus, species)
        
# make function to modify headers
def modify_fna_header(file_name):
    # open contents
    with open(file_name, 'r') as f:
        lines = f.readlines()

    # open the file with writing
    with open(file_name, 'w') as f:
        for line in lines:
            if line.startswith(">"):  # modify only the header lines
                # extract seqid
                seq_id = line.split()[0][1:] # removes first character of line

                if seq_id in seq_info:
                    genus, species = seq_info[seq_id]
                    # modify the header to >SEQ_ID_Genus_Species
                    new_header = f">{seq_id}_{genus}_{species}\n"
                    f.write(new_header)  # make new header
                else:
                    # if not in metadata keep orignal
                    f.write(line)
            else:
                # leave nonheaders unchanged
                f.write(line)

# fun function
for file_name in os.listdir():
    if file_name.endswith(".fna"):
        modify_fna_header(file_name)

print("Headers have been modified successfully.")
