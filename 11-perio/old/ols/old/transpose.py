import pandas as pd
import sys

if len(sys.argv) != 3:
    print("Usage: python script.py inputfile outputfile")
    sys.exit(1)

#import file
input_file = sys.argv[1]  
output_file = sys.argv[2]  

try:
    df = pd.read_csv(input_file, sep='\t')
except FileNotFoundError:
    print(f"Error: The file {input_file} does not exist.")
    sys.exit(1)

#transpose
df_transposed = df.T

df_transposed.reset_index(inplace=True)

df_transposed.to_csv(output_file, sep='\t', index=False)

print(f"Transposed table saved to {output_file}")
