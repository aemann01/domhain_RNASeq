# 1. Get seq ids for ureABC
```sh
cd ~/rna_dohmain/09-urease/01-operon-identify
grep "eC_number=3\.5\.1\.5" ../../homd_map/ALL_genomes.gff| grep ureA | sed 's/|.*//' > ureA.seqs
grep "eC_number=3\.5\.1\.5" ../../homd_map/ALL_genomes.gff| grep ureB | sed 's/|.*//' > ureB.seqs
grep "eC_number=3\.5\.1\.5" ../../homd_map/ALL_genomes.gff| grep ureC | sed 's/|.*//' > ureC.seqs
wc -l *seqs
  # 2489 ureA.seqs
  # 2342 ureB.seqs
  # 2381 ureC.seqs
```
# 2. Find common seqids for ureABC
```sh
python3 common_elements.py #output ureABC.seqs
```
# 3. Get sequence ids to start seeing if they are close enough
```sh
grep "eC_number=3\.5\.1\.5" ../../homd_map/ALL_genomes.gff| grep ureA | sed 's/.*ID=//' | sed 's/;eC.*//'> ureA.seqids
grep "eC_number=3\.5\.1\.5" ../../homd_map/ALL_genomes.gff| grep ureB | sed 's/.*ID=//' | sed 's/;eC.*//'> ureB.seqids
grep "eC_number=3\.5\.1\.5" ../../homd_map/ALL_genomes.gff| grep ureC | sed 's/.*ID=//' | sed 's/;eC.*//'> ureC.seqids

parallel -a ureABC.seqs -j 7 -k "grep '{}' ureA.seqids"> ureA.ids
parallel -a ureABC.seqs -j 7 -k "grep '{}' ureB.seqids"> ureB.ids
parallel -a ureABC.seqs -j 7 -k "grep '{}' ureC.seqids"> ureC.ids
parallel -a ureABC.seqs -j 7 -k "grep -c '{}' ureA.ids"> ureA.counts
parallel -a ureABC.seqs -j 7 -k "grep -c '{}' ureB.ids"> ureB.counts
parallel -a ureABC.seqs -j 7 -k "grep -c '{}' ureC.ids"> ureC.counts
paste -d "," ureABC.seqs ureA.counts ureB.counts ureC.counts >ureABC.counts
```
# 4. Find how many copies different sequences have
```sh
parallel -a ureABC.seqs -j 7 -k "grep '{}' ureA.seqids"> ureA.ids
parallel -a ureABC.seqs -j 7 -k "grep '{}' ureB.seqids"> ureB.ids
parallel -a ureABC.seqs -j 7 -k "grep '{}' ureC.seqids"> ureC.ids
parallel -a ureABC.seqs -j 7 -k "grep -c '{}' ureA.ids"> ureA.counts
parallel -a ureABC.seqs -j 7 -k "grep -c '{}' ureB.ids"> ureB.counts
parallel -a ureABC.seqs -j 7 -k "grep -c '{}' ureC.ids"> ureC.counts
paste -d "," ureABC.seqs ureA.counts ureB.counts ureC.counts >ureABC.counts
```
# 5. Get coordinates
```sh
grep "eC_number=3\.5\.1\.5" ../../homd_map/ALL_genomes.gff| grep ureA > ureA.gff
parallel -a ureA.ids -j 7 -k  "grep '{}' ureA.gff | sed 's/.*CDS\t//'" > ureA.cords
awk -F "\t" '{print $1,$2,$4}' ureA.cords | sed 's/ /\t/g'> ureA.cords2

grep "eC_number=3\.5\.1\.5" ../../homd_map/ALL_genomes.gff| grep ureB > ureB.gff
parallel -a ureB.ids -j 7 -k  "grep '{}' ureB.gff | sed 's/.*CDS\t//'" > ureB.cords
awk '{print $1,$2,$4}' ureB.cords | sed 's/ /\t/g' > ureB.cords2

grep "eC_number=3\.5\.1\.5" ../../homd_map/ALL_genomes.gff| grep ureC > ureC.gff
parallel -a ureC.ids -j 7 -k  "grep '{}' ureC.gff | sed 's/.*CDS\t//'" > ureC.cords
awk '{print $1,$2,$4}' ureC.cords | sed 's/ /\t/g' > ureC.cords2

#clean up files
paste ureA.ids ureA.ids | sed 's/\..*\t/\t/' > tempA
paste -d "\t" tempA ureA.cords2 > ureA.cords

paste ureB.ids ureB.ids | sed 's/\..*\t/\t/' > tempB
paste -d "\t" tempB ureB.cords2 > ureB.cords

paste ureC.ids ureC.ids | sed 's/\..*\t/\t/' > tempC
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
grep "\-40"  ureAB.pos.seqsid | sed 's/\t.*//' | while read line; do grep -m 1 $line ../../homd_map/annotations.merge.txt; done
#the -3s belong to B. subtilis, A. oris, A. naeslundii, and M. tuberculosis
grep "\-3" ureAB.pos.seqsid | sed 's/\t.*//' | while read line; do grep -m 1 $line ../../homd_map/annotations.merge.txt; done

#ureB comes before ureA in the negative strand
awk '$10 ~ /-/ { print $0 }' ureAB.neg.seqsid 
#for some reason all the -40 ones belong to yesinia pestis
grep "\-40"  ureAB.neg.seqsid | sed 's/\t.*//' | while read line; do grep -m 1 $line ../../homd_map/annotations.merge.txt; done
#the -3s belong to B. subtilis, A. oris, A. naeslundii, and M. tuberculosis
grep "\-3" ureAB.neg.seqsid | sed 's/\t.*//' | while read line; do grep -m 1 $line ../../homd_map/annotations.merge.txt; done
#remove the extra column in the negative file
cut --complement -f11 ureAB.neg.seqsid > temp
mv temp ureAB.neg.seqsid
```
# 10. Find ureAB and ureC in an operon
```sh
python3 ureABC_cords.pos.py #out put ureABC.pos.seqsid
#see if there are any negatives
awk '$15 ~ /-10/ { print $0 }' ureABC.pos.seqsid | awk '{print $15}' | sort |uniq 
# -10
# -13
# -3
# -39
# -7
awk '$15 ~ /-39/ { print $0 }' ureABC.pos.seqsid | awk '{print $1}' | while read line; do grep -m 1 $line ../../homd_map/annotations.merge.txt; done
# SEQF3743.1_00001	proB_1	SEQF3743.1	Variovorax	paradoxus
awk '$15 ~ /-7/ { print $0 }' ureABC.pos.seqsid | awk '{print $1}' | while read line; do grep -m 1 $line ../../homd_map/annotations.merge.txt; done
# E. coli, K. pneumonia, E. hormaechei, E. cancerogenus
awk '$15 ~ /-10/ { print $0 }' ureABC.pos.seqsid | awk '{print $1}' | while read line; do grep -m 1 $line ../../homd_map/annotations.merge.txt; done
# SEQF7564.1_00001	none	SEQF7564.1	Staphylococcus	lugdunensis
awk '$15 ~ /-13/ { print $0 }' ureABC.pos.seqsid | awk '{print $1}' | while read line; do grep -m 1 $line ../../homd_map/annotations.merge.txt; done
# Haemophilus influenza (SEQF5362.1) and Corynebacterium bovis
awk '$15 ~ /-3/ { print $0 }' ureABC.pos.seqsid | grep -w "\-3" | awk '{print $1}' | while read line; do grep -m 1 $line ../../homd_map/annotations.merge.txt; done
#lots of different ones

python3 ureABC_cords.neg.py # output ureABC.neg.seqsid
awk '$15 ~ /-10/ { print $0 }' ureABC.neg.seqsid | awk '{print $15}' | sort |uniq 
# -10
awk '$15 ~ /-10/ { print $0 }' ureABC.neg.seqsid | awk '{print $1}' | while read line; do grep -m 1 $line ../../homd_map/annotations.merge.txt; done
# SEQF7533.1_00001	mtrA	SEQF7533.1	Staphylococcus	lugdunensis
# SEQF8590.1_00001	none	SEQF8590.1	Staphylococcus	caprae
# SEQF7554.1_00001	none	SEQF7554.1	Staphylococcus	lugdunensis
#remove the extra column in the negative file
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
# 12. Get some stats for the full operon
```sh
python3 operon.stats
```






