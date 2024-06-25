# 1. Pull sequences from HOMD database
```sh
#remove word wrap from HOMD fna file
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' ../../homd_map/ALL_genomes.fna > ALL_genomes.fna
sed -i 's/|/_/' ALL_genomes.fna
#get specifc sequences
awk '{print $2}' ../01-operon-identify/ureABC.operon | while read line; do grep $line ../01-operon-identify/ureA.gff; done | awk '{print $1, $4, $5}' | sed 's/|/_/' | sort | uniq | sed 's/ /\t/g'> operon.seqs
#getting the genomes
python3 flip_cords2.py #output operon.seq
awk '{print $1}' operon.seq | sort | uniq> operon.seqs2
#pulling the genomes
parallel -a operon.seqs2 -j 7 -k "grep -A 1 '{}' ALL_genomes.fna"> seqs.fa
# 2217 genomes
```
# 2. Get the range of coordinates
```sh
#get back to full ids
sed 's/_.*//' operon.seq | sed 's/\..//' | sed 's/ /\t/g' > temp
awk '{print $1, $2}' operon.seq > operon.seq2
paste -d "\t" temp operon.seq2 | sed 's/ /\t/g' > full.id
python3 combine_seqids.py #output ids.cords and sorted.cords
#get cords in right order and just beginning and end
awk '{print $1, $2, $7}' sorted.cords | sed 's/ /\t/g' | sed '1d'> ids.txt
#subset the operon
# awk -v number="300" '$3!=0{$3+=number} 1' ids.txt > temp #adding 100 the end
# awk -v number="300" '$2!=0{$2+=number} 1' temp > temp1 #adding 100 the end
# mv temp1 ids.txt
seqtk subseq seqs.fa ids.txt > operon.fa
#fix headers
sed 's/ .*//' operon.fa | sed 's/-.*//'> ure_operon.shortid.fasta

#make alignment overall
vsearch --derep_fulllength ure_operon.shortid.fasta \
	--output ure_operon.uniq.fna
mafft --thread 7 \
	--adjustdirectionaccurately \
	ure_operon.uniq.fna > ure_operon.align.fna
	
trimal -in ure_operon.align.fna \
	-out ure_operon.trim.fna \
	-htmlout ure_operon.trim.html \
	-gt 0.5 \
	-resoverlap 0.5 \
	-seqoverlap 50 
#make an alignment for positive ones
grep positive ../01-operon-identify/ureABC.operon | awk '{print $1}' | while read line; do grep $line ids.txt; done > pos.ids
seqtk subseq seqs.fa ids.txt | sed 's/ .*//' operon.fa | sed 's/-.*//' > operon.pos.fa

vsearch --derep_fulllength operon.pos.fa \
	--output ure_operon.pos.uniq.fna
mafft --thread 7 \
	--adjustdirectionaccurately \
	ure_operon.pos.uniq.fna > ure_operon.pos.align.fna

trimal -in ure_operon.pos.align.fna \
	-out ure_operon.pos.trim.fna \
	-htmlout ure_operon.pos.trim.html \
	-gt 0.5 \
	-resoverlap 0.5 \
	-seqoverlap 50 
```
# 3. Make own subsetted gff file
```sh
prokka --cpus 7 \
   --force \
   --prefix ure_operon  \
   --norrna \
   --notrna \
   ure_operon.shortid.fasta \
   1> prokka.out \
   2> prokka.err
#check prokka ouput
gffread -T all_assemblies.gff -o all_assemblies.gtf

grep "product" ure_assemblies.gff | sed 's/.*;pr/pr/' | grep -v Urease | sort | uniq             
# product=Acetyltransferase
# product=hypothetical protein
# product=L-amino acid N-acyltransferase MnaT
# product=L-methionine sulfoximine/L-methionine sulfone acetyltransferase
grep ureC ure_operon.gff | sed 's/:.*//' | while read line; do grep -v $line ure_operon.gff; done | grep SEQ | sed 's/:.*//' | sort | uniq | grep -v "sequence-regi"| grep -v ">" > missing_ureC	
grep ureB ure_operon.gff | sed 's/:.*//' | while read line; do grep -v $line ure_operon.gff; done | grep SEQ | sed 's/:.*//' | sort | uniq | grep -v "sequence-regi"| grep -v ">" > missing_ureB
grep ureA ure_operon.gff | sed 's/:.*//' | while read line; do grep -v $line ure_operon.gff; done | grep SEQ | sed 's/:.*//' | sort | uniq | grep -v "sequence-regi"| grep -v ">" > missing_ureA
```
# 4. Upload filtered files and fna and gtf to palmetto
```sh
rsync -a ~/hiv_rnaseq/merged_Qfiltered scrull@hpcdtn01.rcd.clemson.edu:/scratch/scrull/hiv_rnaseq/filtered
rsync -a ~/rna_dohmain/10-urease/02-operon-mapping/ure_assemblies/ure_assemblies.fna scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/operon-mapping
rsync -a ~/rna_dohmain/10-urease/02-operon-mapping/ure_assemblies/ure_assemblies.gff scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/operon-mapping
```