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
awk -v number="1" '$3!=0{$3+=number} 1' ids.txt > temp #adding 100 the end
mv temp ids.txt
seqtk subseq seqs.fa ids.txt > operon.fa
#fix headers
sed 's/ .*//' operon.fa | sed 's/-.*//'> ure_operon.shortid.fasta

# #make alignment overall
# vsearch --derep_fulllength ure_operon.shortid.fasta \
# 	--output ure_operon.uniq.fna
# mafft --thread 7 \
# 	--adjustdirectionaccurately \
# 	ure_operon.uniq.fna > ure_operon.align.fna
	
# trimal -in ure_operon.align.fna \
# 	-out ure_operon.trim.fna \
# 	-htmlout ure_operon.trim.html \
# 	-gt 0.5 \
# 	-resoverlap 0.5 \
# 	-seqoverlap 50 
# #make an alignment for positive ones
# grep positive ../01-operon-identify/ureABC.operon | awk '{print $1}' | while read line; do grep $line ids.txt; done > pos.ids
# seqtk subseq seqs.fa ids.txt | sed 's/ .*//' operon.fa | sed 's/-.*//' > operon.pos.fa

# vsearch --derep_fulllength operon.pos.fa \
# 	--output ure_operon.pos.uniq.fna
# mafft --thread 7 \
# 	--adjustdirectionaccurately \
# 	ure_operon.pos.uniq.fna > ure_operon.pos.align.fna

# trimal -in ure_operon.pos.align.fna \
# 	-out ure_operon.pos.trim.fna \
# 	-htmlout ure_operon.pos.trim.html \
# 	-gt 0.5 \
# 	-resoverlap 0.5 \
# 	-seqoverlap 50 
```

# 3. Make own subsetted gff file
Run prokka to get template gff file
```sh
prokka --cpus 7 \
   --force \
   --prefix ure_operon  \
   --norrna \
   --notrna \
   ure_operon.shortid.fasta \
   1> prokka.out \
   2> prokka.err
```
```sh
python3 new_cords.py #new.cords
awk '{print $1, $2}' sorted.cords | awk -v number="1" '$2!=0{$2+=number} 1' | sed 's/ /:/' > ids
paste -d "\t" new.cords ids > new.cords2

awk -v OFS='\t' '{print $22, "Prodigal:002006", "CDS", $16, $17, ".", $13, 0, "ID=",$2, ";eC_number=3.5.1.5;Name=ureA;db_xref=COG:COG0831;gene=ureA;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P14916;locus_tag=", $2, ";product=Urease subunit alpha"}' new.cords2 | sed 's/ID=\t/ID=/' | sed 's/tag=\t/tag=/' | sed 's/\t;pro/;pro/' | sed 's/\t;eC/;eC/' | sed 's/positive/+/' | sed 's/negative/-/' | sed '1d'> gff.ureA

awk -v OFS='\t' '{print $22, "Prodigal:002006", "CDS", $18, $19, ".", $13, 0, "ID=",$3, ";eC_number=3.5.1.5;Name=ureB;db_xref=COG:COG0804;gene=ureB;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P69996;locus_tag=", $3, ";product=Urease subunit beta"}' new.cords2 | sed 's/ID=\t/ID=/' | sed 's/tag=\t/tag=/' | sed 's/\t;pro/;pro/' | sed 's/\t;eC/;eC/' | sed 's/positive/+/' | sed 's/negative/-/' | sed '1d' > gff.ureB

awk -v OFS='\t' '{print $22, "Prodigal:002006", "CDS", $20, $21, ".", $13, 0, "ID=",$4, ";eC_number=3.5.1.5;Name=ureC;db_xref=COG:COG0804;gene=ureC;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:Q79VJ3;locus_tag=", $4, ";product=Urease subunit alpha"}' new.cords2 | sed 's/ID=\t/ID=/' | sed 's/tag=\t/tag=/' | sed 's/\t;pro/;pro/' | sed 's/\t;eC/;eC/' | sed 's/positive/+/' | sed 's/negative/-/' | sed '1d'> gff.ureC

grep "sequence-region" ./ure_operon/ure_operon.gff | awk -v number="1" '$4!=0{$4+=number} 1' > header.gff
cat header.gff gff.ureA gff.ureB gff.ureC ure_operon.shortid.fasta > operon.gff
sed -i '1 i\\#\#\gff-version 3' operon.gff
gffread -T operon.gff -o operon.gtf
```
# 4. Map using STAR 
Make star database
```sh
STAR --runMode genomeGenerate \
   --genomeFastaFiles ~/rna_dohmain/10-urease/02-operon-mapping/ure_operon.shortid.fasta \
   --runThreadN 10 \
   --limitGenomeGenerateRAM 228697571594 \
   --sjdbGTFfile ~/rna_dohmain/10-urease/02-operon-mapping/operon.gff \
   --genomeChrBinNbits 10 \
   --limitSjdbInsertNsj 3204994 \
   --genomeSAindexNbases 10\
   --sjdbGTFfeatureExon CDS 2> db.err
```
Make one from prokka output to check
```sh
STAR --runMode genomeGenerate \
   --genomeFastaFiles ~/rna_dohmain/10-urease/02-operon-mapping/ure_operon/ure_operon.fna \
   --runThreadN 10 \
   --limitGenomeGenerateRAM 228697571594 \
   --sjdbGTFfile ~/rna_dohmain/10-urease/02-operon-mapping/ure_operon/ure_operon.gff \
   --genomeChrBinNbits 10 \
   --limitSjdbInsertNsj 3204994 \
   --genomeSAindexNbases 10\
   --sjdbGTFfeatureExon CDS
```
# 5. Upload filtered files and genomeDir to palmetto
```sh
rsync -a ~/hiv_rnaseq/merged_Qfiltered scrull@hpcdtn01.rcd.clemson.edu:/scratch/scrull/hiv_rnaseq/filtered
#GenomeDir
rsync -a ~/rna_dohmain/10-urease/02-operon-mapping/GenomeDir scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/operon-mapping/
rsync -a ~/rna_dohmain/10-urease/02-operon-mapping/operon.gtf scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/operon-mapping/
#prokka
rsync -a ~/rna_dohmain/10-urease/02-operon-mapping/GenomeDir scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/operon-mapping/prokka/
```
# 6. Run STAR mapping on Palmetto
```sh
#!/bin/bash

#SBATCH --job-name operonMap-SAMPLE
#SBATCH --nodes 1
#SBATCH --tasks-per-node 40
#SBATCH --cpus-per-task 1
#SBATCH --mem 750gb
#SBATCH --time 72:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# load module
module add biocontainers star/2.7.10b

# move into scratch
cd /scratch/scrull/hiv_rnaseq/filtered

# run command
STAR --runThreadN 40 \
   --genomeDir /scratch/scrull/hiv_rnaseq/operon-mapping/GenomeDir \
   --readFilesIn <(gunzip -c SAMPLE.1.fastq.gz) <(gunzip -c SAMPLE.2.fastq.gz) \
   --outFileNamePrefix /scratch/scrull/hiv_rnaseq/filtered/SAMPLE.operon \
   --outSAMtype BAM Unsorted \
   --outReadsUnmapped /scratch/scrull/hiv_rnaseq/filtered/SAMPLE.operon.unmapped.bam \
   --quantMode TranscriptomeSAM GeneCounts \
   --alignIntronMax 1 \
   --chimOutType SeparateSAMold
```
Create new sample specific mapping
```sh
ls /scratch/scrull/hiv_rnaseq/filtered | sed 's/\..*//' | sort | uniq | while read line; do sed "s/SAMPLE/$line/g" example.sh > $line.sh ; done
# submit
ls D*sh | while read line; do sbatch $line; done
```
prokka test
```sh
#!/bin/bash

#SBATCH --job-name prokkaMap-SAMPLE
#SBATCH --nodes 1
#SBATCH --tasks-per-node 40
#SBATCH --cpus-per-task 1
#SBATCH --mem 750gb
#SBATCH --time 72:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# load module
module add biocontainers star/2.7.10b

# move into scratch
cd /scratch/scrull/hiv_rnaseq/filtered

# run command
STAR --runThreadN 40 \
   --genomeDir /scratch/scrull/hiv_rnaseq/operon-mapping/prokka/GenomeDir \
   --readFilesIn <(gunzip -c SAMPLE.1.fastq.gz) <(gunzip -c SAMPLE.2.fastq.gz) \
   --outFileNamePrefix /scratch/scrull/hiv_rnaseq/filtered/SAMPLE.prokka \
   --outSAMtype BAM Unsorted \
   --outReadsUnmapped /scratch/scrull/hiv_rnaseq/filtered/SAMPLE.prokka.unmapped.bam \
   --quantMode TranscriptomeSAM GeneCounts \
   --alignIntronMax 1 \
   --chimOutType SeparateSAMold
```
```sh
ls /scratch/scrull/hiv_rnaseq/filtered | sed 's/\..*//' | sort | uniq | while read line; do sed "s/SAMPLE/$line/g" example2.sh > $line.prokka.sh ; done
```
# 7. Build a tree of operon
Build a tree
```sh
#!/bin/bash

#SBATCH --job-name operontree
#SBATCH --nodes 1
#SBATCH --tasks-per-node 60
#SBATCH --cpus-per-task 1
#SBATCH --mem 750gb
#SBATCH --time 72:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# load module
module add mafft/7.471 raxml/8.2.12

# move into scratch
cd /scratch/scrull/hiv_rnaseq/operon-mapping
sed 's/ .*//' operon.fa | sed 's/-.*//' > temp
mv temp operon.fa
sed -i 's/:/-/' operon.fa

mafft --thread 60 \
	--adjustdirectionaccurately \
	operon.fa > operon.align.fa

raxmlHPC-PTHREADS-SSE3 -T 60 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n tre -s operon.align.fa
```
Make tree annotations
```sh
#make annotations file
grep ">" operon.align.fa | sed 's/.*SEQ/SEQ/g' | sed 's/_.*//' | while read line; do grep -wm1 $line ../../homd_map/annotations.merge.txt; done | awk '{print $4, $5}' | sed 's/ /_/' > tre.annots
grep ">" operon.align.fa | sed 's/>//' > tres.seq
paste -d "\t" tres.seq tre.annots > operon_tre.annots
sed -i '1 i\SEQ\ttaxa' operon_tre.annots
#download and look at tree
```
# 8. Make feature counts
```sh
#!/bin/bash

#SBATCH --job-name operonFeature
#SBATCH --nodes 1
#SBATCH --tasks-per-node 60
#SBATCH --cpus-per-task 1
#SBATCH --mem 500gb
#SBATCH --time 24:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

module add subread/2.0.3
cd /scratch/scrull/hiv_rnaseq/filtered
mkdir featurecounts
ls *Aligned.out.bam | sed 's/Aligned.out.bam//' | while read line; do featureCounts -f -p -C -B -M -a ../operon-mapping/operon.gtf -o featurecounts/$line.out -T 60 $line\Aligned.out.bam -t CDS -g transcript_id; done

cd featurecounts
paste *out | grep -v "^#" | awk '{printf "%s\t", $1}{for (i=7;i<=NF;i+=7) printf "%s\t", $i; printf "\n"}' > read_counts.txt
# clean up sample names
sed 's/Aligned.out.bam//g' read_counts.txt | sed 's/.homd//'g > temp
mv temp read_counts.txt
