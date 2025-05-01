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
# move to rpoc folder
cd /home/allie/domhain_RNAseq/11-rpoc_processing
wget https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
unzip new_taxdump.zip
```

### 3. Taxonomic assignment with Kraken2 (HOMD specific database -- made by Suzanne)

```bash
kraken2 --db /home/suzanne/rna_dohmain/rpoc/database/kraken_homd --threads 15 --use-names --output rep_set.kraken.out --unclassified rep_set.unclassified.kraken.out --confidence 0.01 rep_set.fa

# remove unclassified from kraken output to prevent problems downstream
grep -v "^U" rep_set.kraken.out > temp
mv temp rep_set.kraken.out
# get tax ids to query taxonomy from ncbi
awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > query
# get full lineage strings from ncbi
cat query | while read line; do grep -w -m 1 ^$line fullnamelineage.dmp | awk -F"|" '{print $3, $2}' | sed 's/\t//g' | sed 's/  / /g' | sed 's/cellular organisms; //' | sed 's/; /;/g' | sed 's/ /_/g'; done > lineages

awk '{print $2}' rep_set.kraken.out > asvids
paste asvids lineages > taxonomy.txt

grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt

# mafft --thread 55 rep_set.filt.fa > rep_set.align.fa
# raxmlHPC-PTHREADS-SSE3 -T 25 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n ref.tre -s rep_set.align.fa
```

### 4. KMA alternative taxonomic assignment 
NOTE: This does recover some groups (e.g., A. jonsonii) but doesn't have substantially better results when compared to the kraken2 data -- keeping this for legacy, but this is not the data that was used in the paper

Want to test if the loss of some taxa (e.g., Actinomyces sp. Kingella oralis) is due to the kraken2 pipeline. Will be testing an alternative kmer-based alignment program, KMA

Install KMA

```bash
cd
git clone https://bitbucket.org/genomicepidemiology/kma.git
cd kma && make
```

Create database from HOMD (same as used for RNASeq mapping -- STAR) and assign taxonomy

```bash
cd ~/domhain_RNAseq/11-rpoc_processing
~/kma/kma index -i ~/domhain_RNAseq/03-star_map/homd_db/ALL_genomes.fna -o templates !&
~/kma/kma -i rep_set.fa -o kma_out -t_db templates
# this seems to recover some of our missing taxa! Let's see what we can do with it
```

Need to pull taxids from the genbank accessions in the file to get lineage information

```bash
conda activate 2024-HIV_RNASeq
zcat kma_out.frag.gz | awk -F"\t" '{print $6 "\t" $7}' | sed 's/|/\t/' | sed 's/ /\t/' | sed 's/ /_/g' > kma_out.res.query

awk -F"\t" '{print $2}' kma_out.res.query | while read line; do 
    esummary -db nuccore -id "$line" | xtract -pattern DocumentSummary -element Caption,TaxId
    sleep 0.5 # trying to prevent too many requests error
done > taxids.kma &!
# make sure that the values match in each file 
diff <(cut -f2 kma_out.res.query | sed 's/\..*//') <(cut -f1 taxids.kma) # this should come back as nothing

# get lineage information for each NCBI tax ID
awk '{print $2}' taxids.kma | while read line; do grep -w -m 1 ^$line ~/fullnamelineage.dmp | awk -F"|" '{print $3, $2}' | sed 's/\t//g' | sed 's/  / /g' | sed 's/cellular organisms; //' | sed 's/; /;/g' | sed 's/ /_/g'; done > lineages
# merge all results
paste kma_out.res.query taxids.kma lineages > full_kma.res
# pull taxonomy file from full results
awk -F"\t" '{print $4 "\t" $7}' full_kma.res > taxonomy.txt

# filter rep set and clean up
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt
```


















