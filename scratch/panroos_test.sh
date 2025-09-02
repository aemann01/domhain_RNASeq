# cluster using panaroo

panaroo -i ./gff_files/*gff -o results --clean-mode strict -t 190 -a core --aligner mafft -c 0.6 --core_threshold 1.00
numb_refs=$(ls gff_files/*gff | wc -l)

mkdir core_genes
grep -c ">" ./results/aligned_gene_sequences/*fas | grep ":$numb_refs" | sed 's/.\/results\/aligned_gene_sequences\///' | sed 's/:.*//' | while read line; do cp ./results/aligned_gene_sequences/$line ./core_genes/$line; done
cd core_genes
ls ./*fas | sed 's/.\/results\/aligned_gene_sequences\///' | while read line; do grep -n ">" ./$line | sed 's/;.*//' |sort |uniq | wc -l; done | grep "$numb_refs" -c
ls *fas | wc -l
# test for recombination
parallel -j 190 'Phi -o -f {} > {.}.rec' ::: *aln.fas

grep "Normal" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.05)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb1

grep "Max" *.rec | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.05)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb2

grep "NSS" *.rec | grep permutations | while IFS=: read -r file var pval; do
    # Remove leading/trailing whitespace
    pval=$(echo "$pval" | xargs)

    # Check for missing values
    if [[ "$pval" == "--" ]]; then
        echo "$file:$var:        $pval  # Not significant (missing)"
    else
        # Convert scientific notation to float and compare with 0.05
        if awk "BEGIN {exit !($pval < 0.05)}"; then
            echo "$file:$var:        $pval  # Not significant (p > 0.05)"
        fi
    fi
done > recomb3

cat <(awk '{print $1}' recomb1) <(awk '{print $1}' recomb2) <(awk '{print $1}' recomb3)| sed 's/rec.*/fas/' | sort | uniq -c | grep -w 3 | awk '{print $2}' | while read line; do mv $line $line.rd; done

sed -i 's/;.*//' *fas