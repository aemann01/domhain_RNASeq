from Bio import SeqIO
from collections import defaultdict

# in and out files
input_file = "core_genome.align.fa"  
output_file = "combined_core.align.fna"  
# dictionary to store sequences grouped by NZ number
grouped_sequences = defaultdict(str)

# parse the input file
with open(input_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        # extract the NZ number from the header
        header = record.description
        if "NZ_" in header:
            prefix = "NZ_"
        elif "NC_" in header:
            prefix = "NC_"
        else:
            prefix = "UNKNOWN_" 
        nz_number = prefix + header.split("_")[1]  
        sequence = str(record.seq)  # get the sequence as a string
        # combine sequences for the same NZ number
        grouped_sequences[nz_number] += sequence

# write out put file
with open(output_file, "w") as out_handle:
    for nz_number, combined_sequence in grouped_sequences.items():
        new_record = f">{nz_number}\n{combined_sequence}\n"
        out_handle.write(new_record)

#make sure everything works
print(f"Combined sequences have been written to {output_file}")
