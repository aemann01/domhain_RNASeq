# rpoC taxonomic assignment post DADA2 processing

### 1. Install kraken2

```bash
cd # installing in home directory
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh /usr/bin # change to install folder of choice

```

### 2. Install the NCBI taxonomy for cross reference

```bash
wget https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
unzip new_taxdump.zip
# move to rpoc folder
cd /home/allie/domhain_RNAseq/11-rpoc_processing
```


### 2. Taxonomic assignment with Kraken2 (HOMD specific database -- made by Suzanne)

```bash
kraken2 --db /home/suzanne/rna_dohmain/rpoc/database/kraken_homd --threads 15 --use-names --output rep_set.kraken.out --unclassified rep_set.unclassified.kraken.out --confidence 0.01 rep_set.fa

awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > query

cat query | while read line; do grep -w -m 1 ^$line ~/fullnamelineage.dmp | awk -F"|" '{print $3, $2}' | sed 's/\t//g' | sed 's/  / /g' | sed 's/cellular organisms; //' | sed 's/; /;/g' | sed 's/ /_/g'; done > lineages

awk '{print $2}' rep_set.kraken.out > asvids
paste asvids lineages > taxonomy.txt

grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt

mafft --thread 55 rep_set.filt.fa > rep_set.align.fa
raxmlHPC-PTHREADS-SSE3 -T 25 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n ref.tre -s rep_set.align.fa
```










