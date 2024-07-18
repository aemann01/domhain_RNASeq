# 1. Get seq ids for ureABC
```sh
cd ~/rna_dohmain/09-urease/07-denovo-operon
grep "eC_number=3\.5\.1\.5" ../../denovo/all_assemblies.gff| grep ureA | sed 's/\t.*//' > ureA.seqs
grep "eC_number=3\.5\.1\.5" ../../denovo/all_assemblies.gff| grep ureB | sed 's/\t.*//' > ureB.seqs
grep "eC_number=3\.5\.1\.5" ../../denovo/all_assemblies.gff| grep ureC | sed 's/\t.*//' > ureC.seqs
wc -l *seqs
  #  397 ureA.seqs
  #  314 ureB.seqs
  #  418 ureC.seqs
  # 1129 total
```
# 1b. Get denovo taxonomy
Make homd kraken2 database
```sh
mkdir ~/rna_dohmain/denovo/database
cd ~/rna_dohmain/denovo/database
wget https://www.homd.org/ftp/genomes/NCBI/current/fna/ALL_genomes.fna
grep ">" ../..//homd_map/ALL_genomes.fna | sed 's/|.*//' | sed 's/>//'> seqIDs
wget https://www.homd.org/ftp/genomes/PROKKA/current/SEQID_info.txt
sed 's/ /\t/g' SEQID_info.txt | sed 's/\t.*https/\thttps/'| sed 's/\t.*GCA_/\tGCA_/' | sed 's/_[^_]*//2g' > seqIDs2GCA #format it so it is seqid and GCA
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank_historical.txt
sed 's/ /_/g' assembly_summary_genbank.txt | awk '{print $1,$6}' | sed '1d' | sed '1d' | sed 's/ /\t/g' > gca2taxid #get just gca and taxid colum
for i in `cat seqIDs`; do grep $i seqIDs2GCA | awk '{print $2}'; done > GCAs #get the GCA for refrence
sort GCAs| uniq > uniq_GCA
parallel -a uniq_GCA -j 7 -k "grep '{}' gca2taxid"> taxid
#check for missing ids
awk '{print $1}' taxid > tax_acc
awk '{print $2}' seqIDs2GCA | sed '1d'> GCA_rpoc 
cat GCA_rpoc tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids
parallel -a missed_ids -j 7 -k "grep '{}' gca2taxid"> taxid2
cat taxid taxid2 > all_taxids
awk '{print $1}' all_taxids > tax_acc
cat GCA_rpoc tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids
#find missing ID GCA
for i in `cat missed_ids`; do grep -m 1 $i assembly_summary_genbank_historical.txt | awk '{print $1, $6}'; done > taxid3 #get the GCA for refrence
cat taxid taxid2 taxid3 > all_taxids
awk '{print $1}' all_taxids > tax_acc
cat GCA_rpoc tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids
#manually involving stuff
#BE CAREFUL HERE
sed 's/\..[1-2]*//' missed_ids| while read line; do grep -m 1 $line assembly_summary_genbank.txt | awk '{print $1, $6}'; done > taxid4 #get the GCA for refrence
sed  -i '1d' missed_ids
paste missed_ids taxid4 > temp
awk '{print $1, $3}' temp > taxid4
cat taxid taxid2 taxid3 taxid4 | sed 's/ /\t/' > all_taxids
awk '{print $1}' all_taxids > tax_acc
cat GCA_rpoc tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids
#get just uniq GCA
sort all_taxids | uniq > GCA_2_taxid
nano GCA_2_taxid #manually add in missed one 410659
sort seqIDs | uniq > seqs
sort seqs | uniq | wc -l #sanity check
parallel -a seqs -j 7 -k "grep '{}' seqIDs2GCA"> seqf2GCA_seqs
python3 GCA2taxid.py
#get sequences
gffread -x - -g ../../homd_map/ALL_genomes.fna ../../homd_map/ALL_genomes.gff | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" }END { printf "%s", n }' > all_genomes.fasta 

grep ">" all_genomes.fasta | sed 's/>//' > locs
paste <(sed 's/_.*//' locs) locs > locus_tags
python3 locus2taxid.py #output prep_headers.tsv
cat locs <(awk '{print $4}' prep_headers.tsv) | sort | uniq -u | sed 's/_.*//' | sort | uniq # check for any missing
grep ">" all_genomes.fasta -c
wc -l prep_headers.tsv # check matches previous number
awk '{print $3}' prep_headers.tsv > taxids
awk '{print $4}' prep_headers.tsv > locus_tags
#fix the headers
paste locus_tags taxids | sed 's/\t/|kraken:taxid|/' | sed 's/\t/ /' | sed 's/^/>/' > fixed_headers
grep -v ">" all_genomes.fasta  > seqs

paste fixed_headers seqs | sed 's/\t/\n/' > ref.fa

mkdir kraken_homd
kraken2-build --download-taxonomy --db kraken_homd/
kraken2-build --add-to-library ref.fa --db kraken_homd/
kraken2-build --build --max-db-size 8000000000 --db kraken_homd/
```
Assign taxonomy
```sh
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar -xvzf new_taxdump.tar.gz
cd ../
#taxonomic assignment
kraken2 --db ./database/kraken_homd \
	--threads 7 \
	--use-names \
	--output rep_set.kraken.out all_assemblies.fna \
	--unclassified-out rep_set.unclassified.out --confidence 0.01
# 	27332496 sequences (12723.16 Mbp) processed in 318.795s (5144.2 Kseq/m, 2394.61 Mbp/m).
  # 19685599 sequences classified (72.02%)
  # 7646897 sequences unclassified (27.98%
```
nt assign on palmetto
```sh
#!/bin/bash

#SBATCH --job-name denovo_taxonomy
#SBATCH --nodes 1
#SBATCH --tasks-per-node 3
#SBATCH --cpus-per-task 1
#SBATCH --mem 950gb
#SBATCH --time 48:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
cd /scratch/scrull/hiv_rnaseq/denovo-taxonomy
module add kraken2/2.1.2
kraken2 --db ./nt_kraken \
  --threads 7 \
  --use-names \
  --output rep_set.kraken.out all_assemblies.fna \
  --unclassified-out rep_set.unclassified.out --confidence 0.01
#   Loading database information... done.
# 27332496 sequences (12723.16 Mbp) processed in 972.173s (1686.9 Kseq/m, 785.24 Mbp/m).
#   21885818 sequences classified (80.07%)
#   5446678 sequences unclassified (19.93%)
```
USe the output from nt database
```sh
rsync -a scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/denovo-taxonomy/rep_set.kraken.out ./
```
Get full taxonomy
```sh
awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > taxids
sed "s/\t//g" ./database/rankedlineage.dmp > rankedlineage.dmp
sort -t "|" -k 1b,1 rankedlineage.dmp > rankedlineage_sorted
cat rankedlineage_sorted | sed 's/|\{2,\}/|/g' > rankedlineage_clean
# add unclassified to taxonomy file
sed -i '1 i\0|unclassified|' rankedlineage_clean
sed 's/|/\t/' rankedlineage_clean | sed 's/ /_/g' >rankedlineage_clean2
python3 lineage.py #output is lineage
awk -F"\t" '{print $2}' lineage | awk -F\| '{s=$NF;for(i=NF-1;i>=1;i--)s=s FS $i;print s}' | sed 's/^|//' | sed 's/ /_/g' | sed 's/|/;/g' > taxonomy
# merge assemebleies ids and taxonomy
awk '{print $2}' rep_set.kraken.out > denovoids
paste denovoids taxonomy > taxonomy.txt
# fix taxonomy 
python fix_taxonomy.py taxonomy.txt > temp
mv temp taxonomy.txt
sed -i 's/unclassified;Bacteria_unknown;Bacteria_unknown;Bacteria_unknown;Bacteria_unknown;Bacteria_unknown;Bacteria_unknown;Bacteria_unknown/unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified/' taxonomy.txt
sed -i 's/;/\t/g' taxonomy.txt
```
Make annotations file with taxonomy
```sh
awk '{print $1, $2, $3}' all_assemblies.tsv | sed 's/ /_/g' | sed 's/_CDS/\tCDS/' | sed 's/CDS_/CDS\t/' > annotations.denovo.txt 
sed -i 's/locus_tag_ftype_gene/locus_tag\tftype\tgene/' annotations.denovo.txt
#get locus tag to seqid
grep "Name" ./all_assemblies.gff| sed 's/\t.*//' > all.seqs
grep "Name" ./all_assemblies.gff | sed 's/.*ID=//' | sed 's/;eC.*//' | sed 's/;Name.*//'> all.seqids
paste -d "\t" all.seqs all.seqids > all.ids

#combine
python3 merge_annots.py #output annotations.denovo.txt

```
# 2. Find common seqids for ureABC

```sh
python3 common_elements.py #out put is ureABC.seqs
```
# 3. Get sequence ids to start seeing if they are close enough
```sh
grep "eC_number=3\.5\.1\.5" ../../denovo/all_assemblies.gff | grep ureA | sed 's/.*ID=//' | sed 's/;eC.*//'> ureA.seqids
grep "eC_number=3\.5\.1\.5" ../../denovo/all_assemblies.gff | grep ureB | sed 's/.*ID=//' | sed 's/;eC.*//'> ureB.seqids
grep "eC_number=3\.5\.1\.5" ../../denovo/all_assemblies.gff | grep ureC | sed 's/.*ID=//' | sed 's/;eC.*//'> ureC.seqids

paste -d "_" ureA.seqs ureA.seqids > ureA.ids
paste -d "_" ureB.seqs ureB.seqids > ureB.ids
paste -d "_" ureC.seqs ureC.seqids > ureC.ids
```
# 4. Find how many copies different sequences have
```sh
parallel -a ureABC.seqs -j 7 -k "grep -c '{}' ureA.ids"> ureA.counts
parallel -a ureABC.seqs -j 7 -k "grep -c '{}' ureB.ids"> ureB.counts
parallel -a ureABC.seqs -j 7 -k "grep -c '{}' ureC.ids"> ureC.counts
paste -d "," ureABC.seqs ureA.counts ureB.counts ureC.counts >ureABC.counts
```
# 5. Get coordinates
```sh
grep "eC_number=3\.5\.1\.5" ../../denovo/all_assemblies.gff| grep ureA > ureA.gff
parallel -a ureA.seqids -j 7 -k  "grep '{}' ureA.gff | sed 's/.*CDS\t//'" > ureA.cords
awk -F "\t" '{print $1,$2,$4}' ureA.cords | sed 's/ /\t/g'> ureA.cords2

grep "eC_number=3\.5\.1\.5" ../../denovo/all_assemblies.gff| grep ureB > ureB.gff
parallel -a ureB.seqids -j 7 -k  "grep '{}' ureB.gff | sed 's/.*CDS\t//'" > ureB.cords
awk '{print $1,$2,$4}' ureB.cords | sed 's/ /\t/g' > ureB.cords2

grep "eC_number=3\.5\.1\.5" ../../denovo/all_assemblies.gff| grep ureC > ureC.gff
parallel -a ureC.seqids -j 7 -k  "grep '{}' ureC.gff | sed 's/.*CDS\t//'" > ureC.cords
awk '{print $1,$2,$4}' ureC.cords | sed 's/ /\t/g' > ureC.cords2

#clean up files
paste ureA.seqs ureA.ids | sed 's/\..*\t/\t/' > tempA
paste -d "\t" tempA ureA.cords2 > ureA.cords

paste ureB.seqs ureB.ids | sed 's/\..*\t/\t/' > tempB
paste -d "\t" tempB ureB.cords2 > ureB.cords

paste ureC.seqs ureC.ids | sed 's/\..*\t/\t/' > tempC
paste -d "\t" tempC ureC.cords2 > ureC.cords
```
# 6. Check double ureA if they are prokka stutter
```sh
awk -F "," '$2 ~ /2/ { print $0 }' ureABC.counts | sed 's/,.*//' > ureA.doub
cat ureA.doub | while read line; do grep -m 1 $line ureA.cords; done > ureA.cords1
cat ureA.doub | while read line; do grep -m 2 $line ureA.cords | tail -n 1; done > ureA.cords2
python3 ureA_check.py #output
#remake cords file for ureA
grep True ureA.comb | awk '{print $1,$2,$3,$8, $5}' > ureA.comb.filt
#remove the combined sequences from the cords file
cp ureA.cords ureA.cord3

sedscript=$(mktemp)
for i in $(awk '{print $1}' ureA.comb.filt); do
  echo "/$i/d" >> "$sedscript"
done
sed -f "$sedscript" -i ./ureA.cord3
rm -f "$sedscript"

#make cords filt file
cat ureA.cord3 ureA.comb.filt > ureA.cords.filt
```
# 7. Check double ureB if they are prokka stutter
```sh
awk -F "," '$3 ~ /2/ { print $0 }' ureABC.counts | sed 's/,.*//' > ureB.doub
cat ureB.doub | while read line; do grep -m 1 $line ureB.cords; done > ureB.cords1
cat ureB.doub | while read line; do grep -m 2 $line ureB.cords | tail -n 1; done > ureB.cords2
python3 ureB_check.py #output ureB.comb
#remake cords file for ureB
grep True ureB.comb | awk '{print $1,$2,$3,$8, $5}' > ureB.comb.filt
#remove the combined sequences from the cords file
cp ureB.cords ureB.cord3

sedscript=$(mktemp)
for i in $(awk '{print $1}' ureB.comb.filt); do
  echo "/$i/d" >> "$sedscript"
done
sed -f "$sedscript" -i ./ureB.cord3
rm -f "$sedscript"

#make cords filt file
cat ureB.cord3 ureB.comb.filt > ureB.cords.filt
```
# 8. Check double ureC if they are prokka stutter for double and triple
```sh
#double
awk -F "," '$4 ~ /2/ { print $0 }' ureABC.counts | sed 's/,.*//' > ureC.doub
cat ureC.doub | while read line; do grep -m 1 $line ureC.cords; done > ureC.cords1
cat ureC.doub | while read line; do grep -m 2 $line ureC.cords | tail -n 1; done > ureC.cords2
python3 ureC_check.py #output ureC.comb
#remake cords file for ureC
grep True ureC.comb | awk '{print $1,$2,$3,$8, $5}' > ureC.comb.filt1
#remove the combined sequences from the cords file
cp ureC.cords ureC.cord3

sedscript=$(mktemp)
for i in $(awk '{print $1}' ureC.comb.filt); do
  echo "/$i/d" >> "$sedscript"
done
sed -f "$sedscript" -i ./ureC.cord3
rm -f "$sedscript"

#triple
awk -F "," '$4 ~ /3/ { print $0 }' ureABC.counts | sed 's/,.*//' > ureC.doub
cat ureC.doub | while read line; do grep -m 1 $line ureC.cords; done > ureC.cords1
cat ureC.doub | while read line; do grep -m 2 $line ureC.cords | tail -n 1; done > ureC.cords2
cat ureC.doub | while read line; do grep -m 3 $line ureC.cords | tail -n 1; done > ureC.cords3
python3 ureC_check2.py #output ureC.comb
#remake cords file for ureC
grep True ureC.comb2 | awk '{print $1,$2,$3,$12,$5}' > ureC.comb.filt2
#remove the combined sequences from the cords file
cat ureC.comb.filt1 ureC.comb.filt2 > ureC.comb.filt
cp ureC.cords ureC.cord3

sedscript=$(mktemp)
for i in $(awk '{print $1}' ureC.comb.filt); do
  echo "/$i/d" >> "$sedscript"
done
sed -f "$sedscript" -i ./ureC.cord3
rm -f "$sedscript"

#make cords filt file
cat ureC.cord3 ureC.comb.filt > ureC.cords.filt
sed -i 's/ /\t/g' ureC.cords.filt
```
# 9. Find ureA and ureB in an operon
```sh
sed -i 's/ /\t/g' *cords.filt
python3 ureAB_cords.py #ureAB.pos.seqsid and ureAB.neg.seqsid is output
#see if there are any negatives in the positive sense strand
grep "-" ureAB.pos.seqsid
#for some reason all the -40 ones belong to yesinia pestis
grep "\-40"  ureAB.pos.seqsid | sed 's/\t.*//' | while read line; do grep -m 1 $line ../../denovo/annotations.merge.txt; done
#the -3s belong to B. subtilis, A. oris, A. naeslundii, and M. tuberculosis
grep "\-3" ureAB.pos.seqsid | sed 's/\t.*//' | while read line; do grep -m 1 $line ../../denovo/annotations.merge.txt; done

#ureB comes before ureA in the negative strand
awk '$10 ~ /-/ { print $0 }' ureAB.neg.seqsid 
#for some reason all the -40 ones belong to yesinia pestis
grep "\-40"  ureAB.neg.seqsid | sed 's/\t.*//' | while read line; do grep -m 1 $line ../../denovo/annotations.merge.txt; done
#the -3s belong to B. subtilis, A. oris, A. naeslundii, and M. tuberculosis
grep "\-3" ureAB.neg.seqsid | sed 's/\t.*//' | while read line; do grep -m 1 $line ../../denovo/annotations.merge.txt; done
#remove the extra column in the negative file
cut --complement -f11 ureAB.neg.seqsid > temp
mv temp ureAB.neg.seqsid
```
# 10. Find ureAB and ureC in an operon
```sh
python3 ureABC_cords.pos.py #out put ureABC.pos.seqsid
#see if there are any negatives
awk '$15 ~ /-10/ { print $0 }' ureABC.pos.seqsid | awk '{print $15}' | sort |uniq 
python3 ureABC_cords.neg.py # output ureABC.neg.seqsid
awk '$15 ~ /-10/ { print $0 }' ureABC.neg.seqsid | awk '{print $15}' | sort |uniq 
cut --complement -f15 ureABC.neg.seqsid > temp
mv temp ureABC.neg.seqsid
```
# 11. Make master file for ureABC
```sh
#postive
awk '{print $1,$2,$6,$11,$10,$15,$3,$4,$7,$8,$12,$13,$5}' ureABC.pos.seqsid | sed 's/ /\t/g' | sed 's/btwn/btwnab/' | sed 's/strandA/strand/' | sed 's/\+/positive/g' > ureABC.pos
#negative
awk '{print $1,$2,$6,$11,$10,$15,$3,$4,$7,$8,$12,$13,$5}' ureABC.neg.seqsid | sed 's/ /\t/g' | sed 's/btwn/btwnab/' | sed 's/strandA/strand/' | sed 's/\(.*\)-/\1negative/' > ureABC.neg
sed -i '1d' ureABC.neg
cat ureABC.pos ureABC.neg> ureABC.operon
#look at which ones have two full operons
awk '{print $1}' ureABC.operon | sort | uniq -d | while read line; do grep $line ureABC.pos -c ; done
awk '{print $1}' ureABC.operon | sort | uniq -d | while read line; do grep $line ureABC.neg -c ; done
#check for duplicate ureA: should be 0
awk '{print $2}' ureABC.operon | sort | uniq -d | while read line; do grep $line ureABC.pos -c ; done
awk '{print $2}' ureABC.operon | sort | uniq -d | while read line; do grep $line ureABC.neg -c ; done
#check for duplicate ureB: should be 0
awk '{print $3}' ureABC.operon | sort | uniq -d | while read line; do grep $line ureABC.pos -c ; done
awk '{print $3}' ureABC.operon | sort | uniq -d | while read line; do grep $line ureABC.neg -c ; done
#check for duplicate ureB: should be 0
awk '{print $4}' ureABC.operon | sort | uniq -d | while read line; do grep $line ureABC.pos -c ; done
awk '{print $4}' ureABC.operon | sort | uniq -d | while read line; do grep $line ureABC.neg -c ; done
```

# 12. Get subsetted fasta file with operon
```sh
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' ../../denovo/all_assemblies.fna > all_assemblies.fna
#get specifc sequences
awk '{print $2}' ./ureABC.operon | sed 's/.*_//' | while read line; do grep $line ./ureA.gff; done | awk '{print $1, $4, $5}' | sed 's/|/_/' | sort | uniq | sed 's/ /\t/g'> operon.seqs
#getting the genomes
python3 flip_cords2.py #output operon.seq
awk '{print $1}' operon.seq | sort | uniq> operon.seqs2
#pulling the genomes
parallel -a operon.seqs2 -j 7 -k "grep -wA 1 '{}' all_assemblies.fna"> seqs.fa
# 81 genomes