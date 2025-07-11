# 1. Download databases
```sh
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_bacteria.dat.gz
zcat uniprot_trembl_bacteria.dat.gz > uniprot_trembl_bacteria.dat
wget https://raw.githubusercontent.com/tseemann/prokka/master/bin/prokka-uniprot_to_fasta_db
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_bacteria.dat.gz
gunzip uniprot_sprot_bacteria.dat.gz
chmod 755 prokka-uniprot_to_fasta_db
```
# 2. Make database
```sh
sed -i '1 i\//' uniprot_sprot_bacteria.dat
cat uniprot_trembl_bacteria.dat uniprot_sprot_bacteria.dat > uniprot_bacteria.dat
#make prokka db
nohup ./prokka-uniprot_to_fasta_db --verbose --evidence=1 uniprot_bacteria.dat > uniprot
mv uniprot ~/bin/miniconda3/envs/2024-HIV_RNASeq/db/kingdom/Bacteria/sprot

cp ~/bin/miniconda3/envs/2024-HIV_RNASeq/db/kingdom/Bacteria/sprot ./sprot_og

./prokka-uniprot_to_fasta_db --verbose --evidence=1 uniprot_bacteria.dat > uniprot_db 2>uniport_log
cp uniprot_db ~/bin/miniconda3/envs/2024-HIV_RNASeq/db/kingdom/Bacteria/sprot

```
# 3. Run prokka
```sh
sudo docker pull staphb/prokka:latest
awk '!/^>/ { printf "%s", $0; n = "\n" } 
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' ~/rna_dohmain/homd_map/ALL_genomes.fna  > ./ALL_genomes.fna
sed -i 's/ .*//' ./ALL_genomes.fna
sed -i 's/\|/_/' ./ALL_genomes.fna

sudo docker run --rm \
  -v $(pwd):$(pwd) -w $(pwd) \
  staphb/prokka \
  prokka --cpus 60 \
	   --force \
	   --prefix prokka_trembl  \
	   --norrna \
	   --notrna \
	   --proteins ~/rna_dohmain/11-perio/05-TrEMBL/uniprot_db \
	   ~/rna_dohmain/11-perio/05-TrEMBL/ALL_genomes.fna 
```
# 4. Star mappign to refseq genomes for red complex
Download HOMD genomes and run NCBI pgap annotations on them
```sh
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
cd /home/suzanne/rna_dohmain/11-perio/05-TrEMBL
export PGAP_INPUT_DIR=/home/suzanne/docker/.pgap
# now install pgap
./pgap.py --update
exit
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
cd /home/suzanne/rna_dohmain/11-perio/05-TrEMBL
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
cd /home/suzanne/rna_dohmain/11-perio/05-TrEMBL/
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
cd /home/suzanne/rna_dohmain/11-perio/05-TrEMBL/
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

# cd /home/suzanne/rna_dohmain/11-perio/05-TrEMBL
# export TMPDIR=/home/suzanne/tmp
# PGAP_INPUT_DIR=/home/suzanne/docker/.pgap python3 ./pgap.py -r --ignore-all-errors --cpu 60 -o SEQF10074.1_modified -g SEQF10074.1_modified.fna -s 'Lentilactobacillus buchneri' 1> SEQF10074.1_modified.out 2> SEQF10074.1_modified.err
```
After pgap is finished
```sh
for i in ~/rna_dohmain/11-perio/05-TrEMBL/*_modified/annot.gff
do mv $i ${i/\/annot.gff/_annot.gff}
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
rsync -a ~/rna_dohmain/11-perio/05-TrEMBL/combined.fna scrull@hpcdtn01.rcd.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/
rsync -a ~/rna_dohmain/11-perio/05-TrEMBL/combined.gtf scrull@hpcdtn01.rcd.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/
rsync -a ~/rna_dohmain/11-perio/05-TrEMBL/GenomeDir scrull@hpcdtn01.rcd.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/
```
Run star mapping
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







## This is what I actually ended up doing
```sh
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt

paste -d "/" <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.gff.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > trep_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.gff.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > porph_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.gff.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > tann_query

rm query
cat *query > query
wget -i query &
cat <(cat query | sed 's/.*\///') <(ls GCF*) | sort | uniq -u > missing

paste -d "/" <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > trep_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > porph_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > tann_query

rm query
cat *query > query
wget -i query &
cat <(cat query | sed 's/.*\///') <(ls GCF*.fna.gz) | sort | uniq -u > missing
gzip -d GCF*.gz
ls GCF* | sed 's/\..*//' | sort | uniq -u | while read line; do rm $line*; done

paste -d "/" <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.gbff.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > trep_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.gbff.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > porph_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.gbff.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > tann_query

rm query
cat *query > query
wget -i query &
cat <(cat query | sed 's/.*\///') <(ls GCF*) | sort | uniq -u > missing
gzip -d *gz
```
Combine gff file
```py 
from Bio import SeqIO

def merge_gff_files(file_list, output_file):
    with open(output_file, 'w') as outfile:
        for file in file_list:
            with open(file, 'r') as infile:
                for line in infile:
                    if not line.startswith('#'):  # Skip comment lines
                        outfile.write(line)

# List of GFF files
gff_files = [
    "GCF_000007585.1_ASM758v1_genomic.gff",
    "GCF_002754075.1_ASM275407v1_genomic.gff",
    "GCF_024181485.1_ASM2418148v1_genomic.gff",
    "GCF_000008185.1_ASM818v1_genomic.gff",
    "GCF_002754095.1_ASM275409v1_genomic.gff",
    "GCF_024181505.1_ASM2418150v1_genomic.gff",
    "GCF_000010505.1_ASM1050v1_genomic.gff",
    "GCF_002754115.1_ASM275411v1_genomic.gff",
    "GCF_024181545.1_ASM2418154v1_genomic.gff",
    "GCF_000238215.1_ASM23821v1_genomic.gff",
    "GCF_002754135.1_ASM275413v1_genomic.gff",
    "GCF_024181565.1_ASM2418156v1_genomic.gff",
    "GCF_000270225.1_ASM27022v1_genomic.gff",
    "GCF_002754155.1_ASM275415v1_genomic.gff",
    "GCF_024181605.1_ASM2418160v1_genomic.gff",
    "GCF_001263815.1_ASM126381v1_genomic.gff",
    "GCF_002892555.1_ASM289255v1_genomic.gff",
    "GCF_024181625.1_ASM2418162v1_genomic.gff",
    "GCF_001274615.1_ASM127461v1_genomic.gff",
    "GCF_002892575.1_ASM289257v1_genomic.gff",
    "GCF_024181645.1_ASM2418164v1_genomic.gff",
    "GCF_001314265.1_ASM131426v1_genomic.gff",
    "GCF_002892595.1_ASM289259v1_genomic.gff",
    "GCF_024400535.1_ASM2440053v1_genomic.gff",
    "GCF_001444325.1_ASM144432v1_genomic.gff",
    "GCF_018141745.1_ASM1814174v1_genomic.gff",
    "GCF_024400725.1_ASM2440072v1_genomic.gff",
    "GCF_001547855.1_ASM154785v1_genomic.gff",
    "GCF_018141765.1_ASM1814176v1_genomic.gff",
    "GCF_025905525.1_ASM2590552v1_genomic.gff",
    "GCF_001547875.1_ASM154787v1_genomic.gff",
    "GCF_023822425.1_ASM2382242v1_genomic.gff",
    "GCF_028335085.1_ASM2833508v1_genomic.gff",
    "GCF_002753935.1_ASM275393v1_genomic.gff",
    "GCF_023822445.1_ASM2382244v1_genomic.gff",
    "GCF_028335105.1_ASM2833510v1_genomic.gff",
    "GCF_002753955.1_ASM275395v1_genomic.gff",
    "GCF_023822465.1_ASM2382246v1_genomic.gff",
    "GCF_028335125.1_ASM2833512v1_genomic.gff",
    "GCF_002753975.1_ASM275397v1_genomic.gff",
    "GCF_024181405.1_ASM2418140v1_genomic.gff",
    "GCF_030252365.1_ASM3025236v1_genomic.gff",
    "GCF_002754015.1_ASM275401v1_genomic.gff",
    "GCF_024181425.1_ASM2418142v1_genomic.gff",
    "GCF_030440475.1_ASM3044047v1_genomic.gff",
    "GCF_002754035.1_ASM275403v1_genomic.gff",
    "GCF_024181445.1_ASM2418144v1_genomic.gff",
    "GCF_030440495.1_ASM3044049v1_genomic.gff",
    "GCF_002754055.1_ASM275405v1_genomic.gff",
    "GCF_024181465.1_ASM2418146v1_genomic.gff",
    "GCF_900638305.1_57043_C01_genomic.gff",
    "GCF_018141805.1_ASM1814180v1_genomic.gff",
    "GCF_030144345.1_ASM3014434v1_genomic.gff"
]

output_file = "combined.gff"
merge_gff_files(gff_files, output_file)

```
Make db
```sh
cat GCF*fna | sed 's/ .*//' > combined.fna

gffread combined.gff -Tv -o combined.gtf
STAR --runMode genomeGenerate \
	--genomeFastaFiles combined.fna  \
	--runThreadN 7 \
	--sjdbGTFfile combined.gtf \
	--genomeChrBinNbits 10 \
	--limitGenomeGenerateRAM 1493556046773 \
	--sjdbGTFfeatureExon CDS 2> db.log
```
Move db onto palmetto
```sh
#GenomeDir
rsync -a ~/rna_dohmain/filtered scrull@hpcdtn02.rcd.clemson.edu:/scratch/scrull/hiv_rnaseq/

rsync -a ~/rna_dohmain/11-perio/05-TrEMBL/combined.fna scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/
rsync -a ~/rna_dohmain/11-perio/05-TrEMBL/combined.gtf scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/
#prokka
rsync -a ~/rna_dohmain/11-perio/05-TrEMBL/GenomeDir scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/
```
Run star mapping
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
scp  scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/red-mapping/featurecounts/read_counts.txt ~/rna_dohmain/11-perio/05-TrEMBL
```
Make gene map key
```sh
cd ~/rna_dohmain/11-perio/05-TrEMBL
parallel -a <(awk '{print $1}' read_counts.txt | sed '1d') -j 7 -k " grep CDS ./combined.gff | grep -wm 1 '{}'" > gene.txt
#get the contig to taxa
grep CDS ./combined.gff > gene_cds.txt
grep gene= gene_cds.txt > gene2.txt
paste -d "\t" <(awk '{print $1}' gene2.txt) <(sed 's/.*Parent=//' gene2.txt | sed 's/;.*gene=/\t/' | sed 's/;.*product=/\t/' | sed 's/;.*//' | sed 's/ /_/g') > annots.txt
paste -d "/" <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Treponema | grep denticola | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > trep_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Porphyromonas | grep gingivalis | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > porph_query
paste -d "/" <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/') <(cat assembly_summary_refseq.txt | grep Tannerella | grep forsythia | grep "Complete Genome" | awk -F"\t" '{print $20}' | sed 's/\/GFA_/\t\/GFA_/' | awk '{print $0,$NF}' | sed 's/$/_genomic.fna.gz/' | awk '{print $2}' | sed 's/.*GCF/GCF/') > tann_query

rm query
cat *query > query
wget -i query &
cat <(cat query | sed 's/.*\///') <(ls GCF*.fna.gz) | sort | uniq -u > missing
gzip -d GCF*
#get taxonomy
parallel -a <(awk '{print $1}' annots.txt) -j 7 -k '
    grep -wm 1 "{}" ./*genomic.fna > /dev/null && \
    grep -wm 1 "{}" ./*genomic.fna || \
    echo "{} not found"
' > taxa.txt
sed -i 's/not found/Treponema denticola ATCC 35405, complete genome/' taxa.txt

paste -d "\t" annots.txt <(awk -F " " '{print $2}' taxa.txt) <(awk -F " " '{print $2, $3}' taxa.txt | sed 's/ /_/') <(awk -F " " '{print $2, $3, $4, $5}' taxa.txt | sed 's/ /_/g' | sed 's/,_.*//') > red_annots.txt
grep -v cangingivalis red_annots.txt > temp
mv temp red_annots.txt
sed -i '1i\contig\ttag\tgene\tproduct\tgenus\tspecies\tstrain' red_annots.txt
cd ../06-red-complex
#give uniq identifier to read counts
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
' ../05-TrEMBL/read_counts.txt | sed 's/Geneid_1/Geneid/' > read_counts.txt

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
' ../05-TrEMBL/red_annots.txt > red_annots.txt

parallel -a <(awk '{print $2}' ./red_annots.txt | sed 's/tag/Geneid/') -j 7 -k " grep -w '{}' ./read_counts.txt" > red_counts.txt
```
Move onto 06-red-complex
