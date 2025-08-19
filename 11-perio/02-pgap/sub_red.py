import re
import concurrent.futures
import os
from collections import deque

def process_gene_ids(gene_ids, read_counts_file, output_file):
    # Create a deque to hold the results
    output = deque()
    
    # Function to search for gene ID in read count file and collect matches
    def search_for_gene(gene_id):
        pattern = re.compile(rf"\b{re.escape(gene_id)}\b")
        with open(read_counts_file, 'r') as file:
            for line in file:
                if pattern.search(line):
                    output.append(line)
                    break  # Stop once the gene ID is found

    # Use ThreadPoolExecutor for parallel processing
    with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = [executor.submit(search_for_gene, gene_id) for gene_id in gene_ids]

        # Display progress while processing
        total_gene_ids = len(gene_ids)
        processed_count = 0

        for future in concurrent.futures.as_completed(futures):
            processed_count += 1
            if processed_count % 1000 == 0 or processed_count == total_gene_ids:
                print(f"Progress: {processed_count} of {total_gene_ids} gene IDs processed.")

    # Write the results to the output file, optionally sorting them
    with open(output_file, 'w') as output_file_handle:
        output_file_handle.writelines(sorted(output))  # Sorting can be skipped if not needed

    print(f"Processing complete. Results written to {output_file}")

# Example usage
if __name__ == '__main__':
    gene_ids_file = 'red_annots.txt'
    read_counts_file = 'read_counts.txt'
    output_file = 'red_counts.txt'

    # Read gene IDs from the file
    with open(gene_ids_file, 'r') as f:
        gene_ids = [line.split()[1].replace("tag/", "Geneid/") for line in f]

    # Process the gene IDs and read counts
    process_gene_ids(gene_ids, read_counts_file, output_file)
