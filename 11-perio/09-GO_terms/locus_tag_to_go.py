import re

def extract_locus_go_process(file_path, output_file):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    records = []
    block = []
    in_cds = False

    for line in lines:
        # Detect start of a CDS or gene block
        if line.strip().startswith("CDS"):
            in_cds = True
            block = [line]
        elif in_cds and (line.startswith("     ") or line.startswith("                     ")):
            block.append(line)
        else:
            if in_cds:
                # Process the completed block
                locus_tag = None
                go_processes = []

                for bline in block:
                    # Extract locus_tag
                    if "/locus_tag=" in bline:
                        match = re.search(r'/locus_tag="([^"]+)"', bline)
                        if match:
                            locus_tag = match.group(1)

                    # Extract GO_process terms
                    go_matches = re.findall(r'GO_process:\s*(GO:\d+)\s*-\s*([^;\[]+)', bline)
                    for go_id, description in go_matches:
                        go_processes.append(f"{go_id} - {description.strip()}")

                if locus_tag:
                    # If no GO_process found, print "no_term"
                    if not go_processes:
                        records.append((locus_tag, "no_term"))
                    else:
                        for go in go_processes:
                            records.append((locus_tag, go))

                in_cds = False
                block = []

    # Write results to a file in a tabular format
    with open(output_file, 'w') as out_f:
        out_f.write(f"Locus_Tag\tGO_process\n")  # Header
        for locus_tag, go_process in records:
            out_f.write(f"{locus_tag}\t{go_process}\n")

    print(f"Results written to {output_file}")

# Example usage
if __name__ == "__main__":
    input_gbk_file = "combined_homd.gbk"  # Replace with the actual file path if needed
    output_file = "locus2go.txt"
    extract_locus_go_process(input_gbk_file, output_file)
