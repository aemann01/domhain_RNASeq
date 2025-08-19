import pandas as pd

chunk_size = 50000  # Adjust the chunk size based on available memory

# Read files in chunks to avoid memory issues
try:
    # Read large files in chunks
    read_counts_chunks = pd.read_csv("~/rna_dohmain/11-perio/02-pgap/gene_counts.txt", delimiter="\t", chunksize=chunk_size)
    unprot = pd.read_csv("~/rna_dohmain/11-perio/09-GO_terms/locus2go.txt", delimiter="\t")
    goterms = pd.read_csv("~/rna_dohmain/11-perio/09-GO_terms/gene2go.sub", delimiter="\t")

    # Processing chunks of the large file
    for chunk in read_counts_chunks:
        # Process chunk by chunk, merging with smaller dataframes
        chunk['Geneid_base'] = chunk['Geneid'].apply(lambda x: '_'.join(x.split('_')[:2]))
        
        # Merge with unprot dataframe
        merged_chunk = pd.merge(unprot, chunk, right_on="Geneid_base", left_on="Locus_Tag", how="inner")
        
        # You can append to a file as you go or store it in a list to write later
        with open("go_reads.txt", "a") as outfile:
            merged_chunk.to_csv(outfile, sep="\t", index=False, header=outfile.tell() == 0)

    print("Merge completed successfully.")
except Exception as e:
    print(f"Error: {e}")
