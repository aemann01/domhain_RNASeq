#!/usr/bin/python3

from Bio import SeqIO
import pandas as pd

genbank_file = 'combined_genbank_file.gbk'

# need an empty list
data = []

# getting go function that maps to lccus tag
for record in SeqIO.parse(genbank_file, 'genbank'):
    for feature in record.features:
        if feature.type == 'CDS':
            locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
            go_term = feature.qualifiers.get('GO_function', [''])[0]
            data.append({
                'Locus_tag': locus_tag,
                'Go_term': go_term
            })

# make df and save
df = pd.DataFrame(data)
df['Go_term'] = df['Go_term'].str.split(' - ').str[0]
df['Go_term'] = df['Go_term'].replace('', 'no_term')

with open("locus2go.txt", "w") as output:
    df.to_csv(output, sep="\t", index=False)

