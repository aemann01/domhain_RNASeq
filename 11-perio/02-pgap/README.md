# 1. Download HOMD and reannotate using pgap
Download HOMD genomes and run NCBI pgap annotations on them
```sh
cd ~/rna_dohmain/11-perio/02-pgap
# prepare docker
sudo nano /etc/docker/daemon.json 
# # add 
# {
#   "data-root": "/home/suzanne/docker"
# }
sudo rsync -axPS /var/lib/docker/ /home/suzanne/docker
sudo systemctl stop docker
sudo systemctl start docker
sudo docker info | grep 'Docker Root Dir'
sudo docker ps

wget https://raw.githubusercontent.com/ncbi/pgap/refs/heads/master/scripts/pgap.py
sudo chmod +x pgap.py
sudo -i
cd /home/suzanne/rna_dohmain/11-perio/02-pgap
export PGAP_INPUT_DIR=/home/suzanne/docker/.pgap
# now install pgap
./pgap.py --update
exit
mkdir annotations
cd annotations
wget -r -l1 -A 'SEQ*.fna' -nd -np -nH --cut-dirs=3 -q https://www.homd.org/ftp/genomes/PROKKA/current/fna/ | \
  grep -oP "https://www.homd.org/ftp/genomes/PROKKA/current/fna/SEQ.*\.fna" | \
  parallel -j 190 wget -q
wget https://www.homd.org/ftp/genomes/PROKKA/current/SEQID_info.txt
wget https://www.homd.org/ftp/genomes/PROKKA/current/GCA_ID_info.csv
parallel -j 190 "sed -i 's/ /_/g' {}" ::: *.fna
parallel -j 190 "sed -i 's/_.*//' {}" ::: *.fna
parallel -j 190 "sed -i 's/|.*//' {}" ::: *.fna

python3 modify_headers.py
# add unique identifer to each line 
ls *.fna | parallel -j 190 'awk "/^>/ {print \$0 \"_\" ++i} !/^>/ {print \$0}" {} > {/.}_modified.fna'
# get list of too big long file names
ls *_modified.fna | parallel -j 190 "grep '>' {} -m1 | awk '{print length, \$0}' | awk '\$1 >= 51'" | sort -n | sed 's/_.*/_modified.fna/' | sed 's/.*>//' > too_long
cat too_long | while read line; do sed -i 's/ /_/g' $line; done
cat too_long | while read line; do sed -i 's/\.[1234]//g' $line; done
cat too_long | parallel -j 190 "grep '>' {} -m1 | awk '{print length, \$0}' | awk '\$1 >= 51'" | sort -n | sed 's/_.*/_modified.fna/' | sed 's/.*>//' > too_long
cat too_long | sed 's/_.*//' | while read line; do sed -i 's/_subsp.*//' $line.[1234]_modified.fna; done
cat too_long | sed 's/_.*//' | parallel -j 190 "grep '>' {}.[1234]_modified.fna -m1 | awk '{print length, \$0}' | awk '\$1 >= 51'" | sort -n | sed 's/_.*/_modified.fna/' | sed 's/.*>//' > too_long
cat too_long | sed 's/_.*//' | while read line; do sed -i 's/_(TM7)_\[G-1\]_.*//' $line.[1234]_modified.fna; done
cat too_long | sed 's/_.*//' | while read line; do sed -i 's/_\[G-[12345789]\]_bac.*//'  $line.[1234]_modified.fna; done
cat too_long | sed 's/_.*//' | parallel -j 190 "grep '>' {}.[1234]_modified.fna -m1 | awk '{print length, \$0}' | awk '\$1 >= 51'" | sort -n | sed 's/_.*/_modified.fna/' | sed 's/.*>//' > too_long
cat too_long | sed 's/_.*//' | while read line; do sed -i 's/_\[G.*//'  $line.[1234]_modified.fna; done
cat too_long | sed 's/_.*//' | while read line; do sed -i 's/_\[F.*//'  $line.[1234]_modified.fna; done
cat too_long | sed 's/_.*//' | parallel -j 190 "grep '>' {}.[1234]_modified.fna -m1 | awk '{print length, \$0}' | awk '\$1 >= 51'" | sort -n 
# make sure headers are under 50 characters 
ls *_modified.fna | parallel -j 190 "grep '>' {} -m1 | awk '{print length, \$0}'" | sort -n | tail -n 50
# check how the species look
ls *_modified.fna | parallel -j 190 'grep -m 1 ">" {} | awk -F "_" "{print \$2, \$3}" | sed "s/ sp.*//" | sed "s/ \[.*//" | sed "s/ (.*//" | sed "s/  phylum.*//"' | sort | uniq
```
example.sh
```sh
cd /home/suzanne/rna_dohmain/11-perio/02-pgap
export TMPDIR=/home/suzanne/tmp
PGAP_INPUT_DIR=/home/suzanne/docker/.pgap python3 ./pgap.py -r --ignore-all-errors --cpu 190 -o OUT -g SAMPLE -s 'SPEC' 1> OUT.out 2> OUT.err
```
```sh
ls *_modified.fna | while read line; do 
  OUT=$(echo $line | sed 's/.fna//')
  SPEC=$(grep -m 1 ">" $line | awk -F "_" '{print $2, $3}' | sed 's/ sp.*//' | sed 's/ \[.*//' | sed 's/ (.*//' | sed 's/  phylum.*//')
  sed "s/SAMPLE/$line/g" example.sh | sed "s/OUT/$OUT/g" | sed "s/SPEC/$SPEC/g" > $OUT.sh
done

# sudo bash ./SEQF10093.1_modified.sh
sudo -i 
cd /home/suzanne/rna_dohmain/11-perio/02-pgap/
ls *_modified.sh | xargs -n 1 -P 20 bash
exit
```
Get the ones that did not complete a move to hillary, while on hillary
```sh
find ./ -maxdepth 2 -type f -name 'annot.gff' | wc -l
grep ERROR *out -c | grep -v ":0" | sed 's/.out:.*//' > errors
wc -l errors
# see if any og the files with errors contains annots.gff
cat errors | while read dir; do
    if find "$dir" -name "annots.gff" > /dev/null 2>&1; then
        echo "$dir contains annots.gff"
    else
        echo "$dir does not contain annots.gff"
    fi
done > gff_contains
grep -cw not gff_contains

# find directories without gff
for d in *_modified/; do [ ! -f "$d/annot.gff" ] && echo "$d"; done > non.gff
cat non.gff | while read line; do sudo rm -rf $line; done
cat non.gff | sed 's/\//.sh/' | while read line; do sed -i 's/-r/-r --no-internet/' $line; done
sudo -i 
cd /home/suzanne/rna_dohmain/11-perio/02-pgap/
cat non.gff | sed 's/\//.sh/' | xargs -n 1 -P 20 bash
exit

# make unknown organims down to genus lelve
cat non.gff | sed 's/\//.out/' | while read line; do grep -c "Unknown organism" $line; done
grep "Unknown organism" -nc *out | grep ":1" | sed 's/.out.*/.sh/' > unknown
cat unknown | while read file; do sed -i "s/\('[A-Z][a-z]*\) [a-z][a-z\-]*\('\)/\1\2/" "$file"; done
cat unknown |sed  's/.sh//' | while read line; do sudo rm -rf $line; done
cat unknown | xargs -n 1 -P 20 bash
for d in *_modified/; do [ ! -f "$d/annot.gff" ] && echo "$d"; done > non.gff
cat non.gff | sed 's/\//.out/' | while read line; do grep -c "Unknown organism" $line; done
grep "Unknown organism" -nc *out | grep ":1" | sed 's/.out.*/.sh/' > unknown
cat unknown | sed 's/sh/fna/' | while read line; do sed -i 's/ /_/g' $line ; done
cat unknown | while read line; do sed -i 's/TM7 phylum/TM7 phylum sp. oral taxon 352/' $line; done
cat unknown | while read line; do sed -i 's/Saccharibacteria/Candidatus Saccharibacteria bacterium/' $line; done
cat unknown | while read line; do sed -i 's/Gracilibacteria/Candidatus Altimarinus/' $line; done
cat unknown | while read line; do sed -i 's/Absconditabacteria/Candidatus Absconditabacteria/' $line; done
cat unknown | while read line; do sed -i 's/Nanosynbacter/Candidatus Nanosynbacter/' $line; done
cat unknown |sed  's/.sh//' | while read line; do sudo rm -rf $line; done
cat unknown | xargs -n 1 -P 20 bash
# finish the other ones
for d in *_modified/; do [ ! -f "$d/annot.gff" ] && echo "$d"; done > non.gff
cat non.gff | sed 's/\//.fna/' | while read line; do sed -i 's/ /_/g' $line ; done
cat non.gff | while read line; do sudo rm -rf $line; done
cat non.gff | sed 's/\//.sh/' | xargs -n 1 -P 20 bash

for d in *_modified/; do [ ! -f "$d/annot.gff" ] && echo "$d"; done > non.gff
grep "Error: can't get biosample" *out | sed 's/:.*//' | sed 's/out/fna/' | parallel -j 190 'awk "/^>/ {print \$0 \"_\" ++i} !/^>/ {print \$0}" {} > {/.}_modified.fna'
for file in *_modified_modified.fna; do
  mv "$file" "${file/_modified_modified.fna/_modified.fna}"
done
cat non.gff | sed 's/\//.sh/' | while read line; do  sed -i 's/Peptostreptococcaceae /Eubacterium saphenum/' $line; done
cat non.gff | sed 's/\//.fna/' | while read line; do sed -i 's/\[//' $line; done 
cat non.gff | sed 's/\//.fna/' | while read line; do sed -i 's/\]//' $line; done 
cat non.gff | while read line; do sudo rm -rf $line; done
cat non.gff | sed 's/\//.sh/' | xargs -n 1 -P 20 bash

grep "no such rank in the" -nc *out | grep -v ":0" | sed 's/.out.*/.sh/' > unknown


#see how many CDS each genome has
find ./ -maxdepth 2 -type f -name 'annot.gff' | while read line; do grep CDS -c $line; done 

grep "no such rank in" -nc *out | grep -v ":0" | sed 's/.out.*/.sh/' > no_rank
cat no_rank | while read line; do sed -i 's/Candidatus Saccharibacteria/Candidatus Saccharibacteria bacterium/' $line; done

# cd /home/suzanne/rna_dohmain/11-perio/02-pgap
# export TMPDIR=/home/suzanne/tmp
# PGAP_INPUT_DIR=/home/suzanne/docker/.pgap python3 ./pgap.py -r --ignore-all-errors --cpu 60 -o SEQF10074.1_modified -g SEQF10074.1_modified.fna -s 'Lentilactobacillus buchneri' 1> SEQF10074.1_modified.out 2> SEQF10074.1_modified.err
```
After pgap is finished
```sh
for i in ~/rna_dohmain/11-perio/02-pgap/annotations/*_modified/annot_cds_from_genomic.fna
do mv $i ${i/\/annot_cds_from_genomic.fna/_annot_cds_from_genomic.fna}
done
for i in ~/rna_dohmain/11-perio/02-pgap/annotations/*_modified/annot.gff
do mv $i ${i/\/annot.gff/_annot.gff}
done
for i in ~/rna_dohmain/11-perio/02-pgap/annotations/*_modified/annot.gbk
do mv $i ${i/\/annot.gbk/_annot.gbk}
done
python3 combine_gff.py 
ls *modified_annot.gff | sed 's/_annot.gff/.fna/' | while read line; do cat $line; done > combined.fna
# make db
mkdir gff2gtf
cp combined.gff gff2gtf
split -l 10000 ./gff2gtf/combined.gff ./gff2gtf/combined_chunk_
ls ./gff2gtf/combined_chunk_* | parallel "gffread {} -Tv -o {.}.gtf"
cat ./gff2gtf/combined_chunk_*.gtf > combined.gtf

STAR --runMode genomeGenerate \
    --genomeFastaFiles combined.fna  \
    --runThreadN 190 \
    --sjdbGTFfile combined.gtf \
    --genomeChrBinNbits 10 \
    --limitGenomeGenerateRAM 1493556046773 \
    --sjdbGTFfeatureExon CDS 2> db.log

#GenomeDir
rsync -a ~/rna_dohmain/filtered scrull@hpcdtn01.rcd.clemson.edu:/scratch/scrull/hiv_rnaseq/
rsync -a ~/rna_dohmain/11-perio/02-pgap/combined.fna scrull@hpcdtn01.rcd.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/
rsync -a ~/rna_dohmain/11-perio/02-pgap/combined.gtf scrull@hpcdtn01.rcd.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/
rsync -a ~/rna_dohmain/11-perio/02-pgap/GenomeDir scrull@hpcdtn01.rcd.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/
```
# 4. Run star mapping: example.sh
```sh
#!/bin/bash

#SBATCH --job-name redmap-SAMPLE
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
   --genomeDir /scratch/scrull/hiv_rnaseq/red-mapping/GenomeDir \
   --readFilesIn <(gunzip -c SAMPLE.1.fastq.gz) <(gunzip -c SAMPLE.2.fastq.gz) \
   --outFileNamePrefix /scratch/scrull/hiv_rnaseq/red-mapping/SAMPLE.red \
   --outSAMtype BAM Unsorted \
   --outReadsUnmapped /scratch/scrull/hiv_rnaseq/red-mapping/SAMPLE.red.unmapped.bam \
   --quantMode TranscriptomeSAM GeneCounts \
   --alignIntronMax 1 \
   --chimOutType SeparateSAMold
```
Creapt scripts and submit them
```sh
ls /scratch/scrull/hiv_rnaseq/filtered | sed 's/\..*//' | sort | uniq | while read line; do sed "s/SAMPLE/$line/g" example.sh > $line.sh ; done

ls /home/scrull/rna_scripts/02-operon-mapping/DM*.sh | sed 's/\..*//' | sort | uniq | sed 's/.*\///' | sort | uniq | while read line; do sed "s/SAMPLE/$line/g" example.sh > $line.sh ; done
# submit
ls D*sh | while read line; do sbatch $line; done
```
Make feauture counts
```sh
#!/bin/bash

#SBATCH --job-name redFeature
#SBATCH --nodes 1
#SBATCH --tasks-per-node 60
#SBATCH --cpus-per-task 1
#SBATCH --mem 500gb
#SBATCH --time 24:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

module add subread/2.0.3
cd /scratch/scrull/hiv_rnaseq/red-mapping
mkdir featurecounts
ls *Aligned.out.bam | sed 's/Aligned.out.bam//' | while read line; do featureCounts -f -p -C -B -M -a ./combined.gtf -o featurecounts/$line.out -T 60 $line\Aligned.out.bam -t CDS -g transcript_id; done

cd featurecounts
paste *out | grep -v "^#" | awk '{printf "%s\t", $1}{for (i=7;i<=NF;i+=7) printf "%s\t", $i; printf "\n"}' > read_counts.txt
# clean up sample names
sed 's/Aligned.out.bam//g' read_counts.txt | sed 's/.homd//'g > temp
mv temp read_counts.txt
```
Move feature counts to stella
```sh
scp scrull@hpcdtn01.rcd.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/featurecounts/read_counts.txt ~/rna_dohmain/11-perio/02-pgap
```
# 5. Make annotation file
```sh
cd ~/rna_dohmain/11-perio/02-pgap
awk '
{
    # Initialize or increment the counter for each unique gene ID
    if (!seen[$1]) {
        seen[$1] = 1
    } else {
        seen[$1]++
    }
    
    # Print the updated line with the new suffix
    printf("%s_%d", $1, seen[$1])
    for (i = 2; i <= NF; i++) {
        printf("\t%s", $i)
    }
    print ""
}
' ./read_counts.txt | sed 's/Geneid_1/Geneid/' > read_counts2.txt
# parallel -a <(awk '{print $1}' read_counts.txt | sed '1d') -j 190 -k " grep CDS ./combined.gff | grep -wm 1 '{}'" > gene.txt
#get the genome to taxa
grep CDS ./combined.gff > gene_cds.txt
grep gene= gene_cds.txt > gene2.txt
paste -d "\t" <(awk '{print $1}' gene2.txt) <(sed 's/.*Parent=//' gene2.txt | sed 's/;.*gene=/\t/' | sed 's/;.*product=/\t/' | sed 's/;.*//' | sed 's/ /_/g') > annots.txt

#get taxonomy
python3 get_taxa.py
sed -i '1d' taxa.csv
paste -d "\t" annots.txt <(awk -F "," '{print $3}' taxa.csv) <(awk -F "," '{print $3, $4}' taxa.csv | sed 's/ /_/') <(awk -F "," '{print $3, $4, $5}' taxa.csv | sed 's/ /_/g' | sed 's/,_.*//') > gene_annots.txt
sed -i '1i\genome\ttag\tgene\tproduct\tgenus\tspecies\tstrain' gene_annots.txt
#gives unique identifered to annotations
awk '
BEGIN { FS = OFS = "\t" }  # Set input and output field separators to tab
NR == 1 { print; next }    # Print the header line unchanged
{
    if (!seen[$2]) {
        seen[$2] = 1
    } else {
        seen[$2]++
    }
    
    # Print the modified line with the incremented suffix in the tag column
    $2 = $2 "_" seen[$2]
    print
}
' ./gene_annots.txt > clean_annots.txt
mv ./gene_annots.txt ./gene_annots.txt2 && mv ./clean_annots.txt ./gene_annots.txt
#get just the counts that have CDS
mono-csc sub_annots.cs -out:sub_annots.exe
mono sub_annots.exe
cat <(head read_counts.txt -n 1) gene_counts.txt > temp
mv temp gene_counts.txt
# parallel -a <(awk '{print $2}' ./gene_annots.txt | sed 's/tag/Geneid/') -j 190 -k " grep -wm 1 '{}' ./read_counts.txt" > gene_counts.txt
```
# 6. Just get red complex
```sh
grep -i "Treponema_denticola\|Porphyromonas_gingivalis\|Tannerella_forsythia" clean_annots.txt > red_annots.txt
sed -i '1i\genome\ttag\tgene\tproduct\tgenus\tspecies\tstrain' red_annots.txt
mono-csc sub_red.cs -out:sub_red.exe
mono sub_red.exe
cat <(head gene_counts.txt -n 1) red_counts.txt > temp
mv temp red_counts.txt
# parallel -a <(awk '{print $2}' ./red_annots.txt | sed 's/tag/Geneid/') -j 190 -k " grep -wm 1 '{}' ./read_counts.txt" > red_counts.txt
# cat <(head -n 1 read_counts.txt) red_counts.txt > temp
# mv temp red_counts.txt
```