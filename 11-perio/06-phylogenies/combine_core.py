import os
from collections import defaultdict

# get the current directory
current_dir = os.getcwd()
alignment_folder = current_dir

# Dictionary to store sequences based on their headers
sequences = defaultdict(str)

# Loop through all .align.fa files in the specified folder
for filename in os.listdir(alignment_folder):
    if filename.endswith(".align.fa"):
        file_path = os.path.join(alignment_folder, filename)
        
        # Open each FASTA file
        with open(file_path, "r") as f:
            current_header = None
            current_sequence = []

            # Loop through each line in the file
            for line in f:
                line = line.strip()
                
                # If the line starts with a header (i.e., '>')
                if line.startswith(">"):
                    if current_header and current_sequence:
                        # Add the sequence to the dictionary under the current header
                        sequences[current_header] += ''.join(current_sequence)

                    # Reset for the next sequence
                    current_header = line[1:]  # Remove the '>' character
                    current_sequence = []
                else:
                    current_sequence.append(line)

            # Don't forget to add the last sequence if the file ends
            if current_header and current_sequence:
                sequences[current_header] += ''.join(current_sequence)

# Write the combined sequences to a new FASTA file
output_file = "core_genome.align.fa"
with open(output_file, "w") as out_f:
    for header, seq in sequences.items():
        out_f.write(f">{header}\n")
        out_f.write(seq + "\n")

print(f"Combined sequences written to {output_file}")
