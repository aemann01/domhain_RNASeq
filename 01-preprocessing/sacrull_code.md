# 1. Make directory with complete files
```sh
qsub -I -l select=1:ncpus=2:mem=50gb:interconnect=1g,walltime=4:00:00
find . -type f -exec touch {} +
sudo rsync -a /home/allie/domhain_RNASeq-1.17.24 scrull@hpcdtn02.rcd.clemson.edu:/scratch/scrull/hiv_rnaseq
for i in `ls *L001*.fastq.gz`; do mv "$i" "`echo $i | tr 'L001' 'L009'`"; done # to prevent same sample names. MAKE SURE SAMPLES DO NOT OVERWRITE


#PBS -N copy_files
#PBS -l select=1:ncpus=10:mem=100gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe
cd /scratch/scrull/hiv_rnaseq/

cp NS3597-VRichards*/*.fastq.gz ./data
```
# 2. Fastqc
```sh
cd /scratch1/scrull/hiv_rnaseq/
mkdir fastqc
cp NS3597-VRichards*/*fastq.gz ./fastqc
fastqc -o fastqc *fastq.gz

#PBS -N fastqc
#PBS -l select=1:ncpus=30:mem=100gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe
source ~/.bashrc
conda activate hiv_rnaseq

#module add multiqc/1.12
cd /scratch1/scrull/hiv_rnaseq/data
mkdir ../fastqc
fastqc -o ../fastqc *fastq.gz
```
# 3. Cutadapt
```sh
#PBS -N cutadapt
#PBS -l select=1:ncpus=10:mem=300gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe

# move to folder with data
# CHANGE TO DATA FOLDER
source ~/.bashrc
cd /scratch1/scrull/hiv_rnaseq/data
mkdir cutadapt

# load modules
conda activate hiv_rnaseq

# run
ls ./*_R1_*.fastq.gz | sed 's/_R1_001.fastq.gz//' | while read line; do
	cutadapt \
	-a AGATCGGAAGAG \
	-A AGATCGGAAGAG \
	--nextseq-trim=20 \
	-o cutadapt/$line.1.trim.fastq.gz \
	-p cutadapt/$line.2.trim.fastq.gz \
	--trim-n \
	--minimum-length 100 \
	--max-n 0 \
	-q 24,24 \
	--cores=10 \
	$line\_R1_001.fastq.gz \
	$line\_R2_001.fastq.gz \
	1>cutadapt/$line.trim.out;
	done
```
# 4. Merge fastq files together
```sh
#get sample list
cd /scratch1/scrull/hiv_rnaseq/data/cutadapt
ls *out | sed 's/_.*//' | sort | uniq > /scratch1/scrull/hiv_rnaseq/data/samp.ids
# pbs script for running on palmetto
#PBS -N merge
#PBS -l select=1:ncpus=3:mem=20gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe

# move to folder with data
cd /scratch1/scrull/hiv_rnaseq/
mkdir merged_Qfiltered

# run
cat ./data/samp.ids | while read line; do cat ./data/cutadapt/$line*.1.trim.fastq.gz > ./merged_Qfiltered/$line.1.fastq.gz; done
cat ./datasamp.ids | while read line; do cat ./data/cutadapt/$line*.2.trim.fastq.gz > ./merged_Qfiltered/$line.2.fastq.gz; done
```
# 5. Check Merge files
Stole from Allie since scratch1 and fastscratch were down
```sh
ls DM00346V2PQ54*1.trim.fastq.gz | while read line; do zcat $i | wc -l; done
```
# 6. Remove human/rRNA contamination
First, download latest release of the silva database, concatentate together
```sh
mkdir /scratch1/scrull/hiv_rnaseq/rRNA_map
cd /scratch1/scrull/hiv_rnaseq/rRNA_map
# download latest version of silva LSU and SSU database
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
gzip -d SILVA_138.1_*
cat SILVA_138.1_LSURef_NR99_tax_silva.fasta SILVA_138.1_SSURef_NR99_tax_silva.fasta > silva_rRNA.fa
rm SILVA_138.1_*
```
Download human genome
```sh
cd ../
mkdir human_map
cd /scratch1/scrull/hiv_rnaseq/human_map
# download human reference genome 
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh38.p14.genome.fa.gz
gzip -d GRCh38.p14.genome.fa.gz
```
Index database build for bowtie (~3 hours)
```sh
#PBS -N bwt_index
#PBS -l select=1:ncpus=2:mem=20gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe

# move to folder with data
source ~/.bashrc
conda activate hiv_rnaseq
cd /scratch1/scrull/hiv_rnaseq/rRNA_map

# run
bowtie2-build silva_rRNA.fa silva_rRNA.db

# do same for human genome
cd /scratch1/scrull/hiv_rnaseq/human_map
bowtie2-build GRCh38.p14.genome.fa GRCh38.p14.genome.db
```
Map data to humnan genome
```sh
#PBS -N human_map
#PBS -l select=1:ncpus=40:mem=1500gb
#PBS -l walltime=138:00:00
#PBS -j oe
#PBS -m abe
#PBS -q bigmem

# move to folder with data
source ~/.bashrc
conda activate hiv_rnaseq
cd /scratch/scrull/hiv_rnaseq/merged_Qfiltered

# run
ls *.1.fastq.gz | sed 's/.1.fastq.gz//' | while read line; do
	bowtie2 -x /scratch1/scrull/hiv_rnaseq/human_map/GRCh38.p14.genome.db \
	-1 $line.1.fastq.gz \
	-2 $line.2.fastq.gz \
	--end-to-end  \
	--qc-filter \
	--no-unal \
	--no-head \
	--no-sq \
	-t \
	--threads 40 \
	-S /scratch1/scrull/hiv_rnaseq/human_map/$line.sam \
	2>/scratch1/scrull/hiv_rnaseq/human_map/$line.out \
	1>/scratch1/scrull/hiv_rnaseq/human_map/$line.err;
	done
```
Map to SILVA
```sh
#PBS -N silva_map
#PBS -l select=1:ncpus=40:mem=1500gb
#PBS -l walltime=138:00:00
#PBS -j oe
#PBS -m abe
#PBS -q bigmem

# move to folder with data
source ~/.bashrc
conda activate hiv_rnaseq
cd /scratch/scrull/hiv_rnaseq/merged_Qfiltered

# run
ls *.1.fastq.gz | sed 's/.1.fastq.gz//' | while read line; do
	bowtie2 -x /scratch1/amann3/domhain_RNASeq-1.17.24/rRNA_map/silva_rRNA.db \
	-1 $line.1.fastq.gz \
	-2 $line.2.fastq.gz \
	--end-to-end  \
	--qc-filter \
	--no-unal \
	--no-head \
	--no-sq \
	-t \
	--threads 40 \
	-S /scratch1/scrull/hiv_rnaseq/rRNA_map/$line.sam \
	2>/scratch1/scrull/hiv_rnaseq/rRNA_map/$line.out \
	1>/scratch1/scrull/hiv_rnaseq/rRNA_map/$line.err;
	done
```
Get sequence ids and remove
```sh
#PBS -N seq_ids
#PBS -l select=1:ncpus=5:mem=700gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe

# move to folder with data
source ~/.bashrc
conda activate hiv_rnaseq
cd /scratch1/scrull/hiv_rnaseq/human_map
ls *sam | while read line; do awk '{print $1}' $line | sort | uniq > $line.ids; done
cd ../rRNA_map
ls *sam | while read line; do awk '{print $1}' $line | sort | uniq > $line.ids; done
```
Cat sequences to remove from each sample
```sh
cd /scratch1/scrull/hiv_rnaseq/human_map
ls *ids | sed 's/.sam.ids//' | parallel 'cat {}.sam.ids ../rRNA_map/{}.sam.ids | sort | uniq > ../filtered/{}.filt.ids'
```
Remove sequences from merged files
```sh
#PBS -N remove_seqs
#PBS -l select=1:ncpus=45:mem=1500gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe
#PBS -q bigmem

# move to folder with data
source ~/.bashrc
conda activate hiv_rnaseq
cd /scratch1/scrull/hiv_rnaseq/filtered
ls *ids | sed 's/.filt.ids//' | parallel --gnu -j 44 'python ~/rna_scripts/remove_seqs.py -f /scratch/scrull/hiv_rnaseq/merged_Qfiltered/{}.1.fastq.gz -i {}.filt.ids -o /scratch1/scrull/hiv_rnaseq/filtered/{}.1.fastq.gz'
ls *ids | sed 's/.filt.ids//' | parallel --gnu -j 44 'python ~/rna_scripts/remove_seqs.py -f /scratch/scrull/hiv_rnaseq/merged_Qfiltered/{}.2.fastq.gz -i {}.filt.ids -o /scratch1/scrull/hiv_rnaseq/filtered/{}.2.fastq.gz'
find . -depth -name "*.fastq.gz" -exec sh -c 'f="{}"; mv -- "$f" "${f%.fastq.gz}.fastq"' \;
gzip *fastq
```
# 6. Denovo assembey
```sh
#PBS -N denovo_DM00008V1PQ16-2
#PBS -l select=1:ncpus=41:mem=1500gb
#PBS -l walltime=168:00:00
#PBS -j oe
#PBS -m abe
#PBS -q bigmem

# move to folder with data
cd /scratch1/scrull/hiv_rnaseq/filtered/
source ~/.bashrc
conda activate hiv_rnaseq
# run
   Trinity \
   --seqType fq \
   --max_memory 750G \
   --left DM00008V1PQ16-2.1.fastq.gz \
   --right DM00008V1PQ16-2.2.fastq.gz \
   --CPU 40 \
   --output /scratch1/scrull/denovo/trinity_denovo_DM00008V1PQ16-2 \
   --full_cleanup
```
example.pbs
```sh
#PBS -N denovo_SAMP
#PBS -l select=1:ncpus=41:mem=760gb
#PBS -l walltime=168:00:00
#PBS -j oe
#PBS -m abe
#PBS -q bigmem

# move to folder with data
cd /scratch1/scrull/hiv_rnaseq/filtered/
source ~/.bashrc
conda activate hiv_rnaseq
# run
   Trinity \
   --seqType fq \
   --max_memory 750G \
   --left SAMP.1.fastq.gz \
   --right SAMP.2.fastq.gz \
   --CPU 40 \
   --no_salmon \
   --output /scratch1/scrull/hiv_rnaseq/denovo/trinity_denovo_SAMP \
   --full_cleanup
```
```sh
ls /scratch1/scrull/hiv_rnaseq/filtered/*fastq* | sed 's/\..*//' | sed 's/.*DM/DM/' | sort | uniq > sample.list
cat sample.list | while read line; do sed "s/SAMP/$line/g" example.pbs > $line.pbs; done
ls *.pbs | while read line; do qsub $line; done
```
