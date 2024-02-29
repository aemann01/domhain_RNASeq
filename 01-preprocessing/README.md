# Read Preprocessing

1. Download raw fastq files from /////, activate environment

```bash
wget //////////
# currently running on palmetto
```

2. Run fastqc to generate quality metrics for each sample

```bash
mkdir fastqc
fastqc -o fastqc *fastq.gz
```

3. Next run cutadapt to remove TruSeq adapters, poly G tails, and filtering by quality score with a minimum of Q30

```bash
#PBS -N cutadapt
#PBS -l select=1:ncpus=10:mem=300gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe

# move to folder with data
# CHANGE TO DATA FOLDER
cd /scratch1/amann3/domhain_RNASeq-1.17.24/NS3597-VRichards_10B-22H5YLLT3-Lane1-8
mkdir cutadapt

# load modules
module add cutadapt/4.1

# run
ls *_R1_*.fastq.gz | sed 's/_R1_001.fastq.gz//' | while read line; do
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

4. Merge fastq files from separate lanes into single R1 and R2 files

```bash
# first get a list of sample names to iterate through
ls *out | sed 's/_.*//' | sort | uniq > /scratch1/amann3/domhain_RNASeq-1.17.24/samp.ids
# pbs script for running on palmetto
#PBS -N merge
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe

# move to folder with data
cd /scratch1/amann3/domhain_RNASeq-1.17.24/
mkdir merged_Qfiltered

# run
cat samp.ids | while read line; do cat NS3597-VRichards*/cutadapt/$line*.1.trim.fastq.gz > merged_Qfiltered/$line.1.fastq.gz; done
cat samp.ids | while read line; do cat NS3597-VRichards*/cutadapt/$line*.2.trim.fastq.gz > merged_Qfiltered/$line.2.fastq.gz; done
```

5. Remove human/rRNA contamination

First, download latest release of the silva database, concatentate together

```bash
cd /scratch1/amann3/domhain_RNASeq-1.17.24/rRNA_map
# download latest version of silva LSU and SSU database
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
gzip -d SILVA_138.1_*
cat SILVA_138.1_LSURef_NR99_tax_silva.fasta SILVA_138.1_SSURef_NR99_tax_silva.fasta > silva_rRNA.fa
rm SILVA_138.1_*
```

Download human reference genome for mappint

```bash
cd /scratch1/amann3/domhain_RNASeq-1.17.24/human_map
# download human reference genome 
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh38.p14.genome.fa.gz
gzip -d GRCh38.p14.genome.fa.gz
```

Index databases for bowtie2

```bash
#PBS -N bwt_index
#PBS -l select=1:ncpus=1:mem=20gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe

# move to folder with data
cd /scratch1/amann3/domhain_RNASeq-1.17.24/rRNA_map

# load modules
module add bowtie2/2.5.2

# run
bowtie2-build silva_rRNA.fa silva_rRNA.db

# do same for human genome
cd /scratch1/amann3/domhain_RNASeq-1.17.24/human_map
bowtie2-build GRCh38.p14.genome.fa GRCh38.p14.genome.db
```

Map merged data to human genome

```bash
#PBS -N human_map
#PBS -l select=1:ncpus=40:mem=750gb
#PBS -l walltime=138:00:00
#PBS -j oe
#PBS -m abe
#PBS -q bigmem

# move to folder with data
cd /scratch1/amann3/domhain_RNASeq-1.17.24/merged_Qfiltered

# load modules
module add bowtie2/2.5.2

# run
ls *.1.fastq.gz | sed 's/.1.fastq.gz//' | while read line; do
	bowtie2 -x /scratch1/amann3/domhain_RNASeq-1.17.24/human_map/GRCh38.p14.genome.db \
	-1 $line.1.fastq.gz \
	-2 $line.2.fastq.gz \
	--end-to-end  \
	--qc-filter \
	--no-unal \
	--no-head \
	--no-sq \
	-t \
	--threads 40 \
	-S /scratch1/amann3/domhain_RNASeq-1.17.24/human_map/$line.sam \
	2>/scratch1/amann3/domhain_RNASeq-1.17.24/human_map/$line.out \
	1>/scratch1/amann3/domhain_RNASeq-1.17.24/human_map/$line.err;
	done
```

Map merged data to SILVA

```bash
#PBS -N silva_map
#PBS -l select=1:ncpus=40:mem=750gb
#PBS -l walltime=138:00:00
#PBS -j oe
#PBS -m abe
#PBS -q bigmem

# move to folder with data
cd /scratch1/amann3/domhain_RNASeq-1.17.24/merged_Qfiltered

# load modules
module add bowtie2/2.5.2

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
	-S /scratch1/amann3/domhain_RNASeq-1.17.24/rRNA_map/$line.sam \
	2>/scratch1/amann3/domhain_RNASeq-1.17.24/rRNA_map/$line.out \
	1>/scratch1/amann3/domhain_RNASeq-1.17.24/rRNA_map/$line.err;
	done
```

8. Remove rRNA and human sequences

First, get list of ids to remove

```bash
#PBS -N filt_ids
#PBS -l select=1:ncpus=10:mem=750gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe
#PBS -q bigmem

# move to folder with data
cd /scratch/amann3/human_map

# get ids that mapped and remove sam files to free up space
ls *sam | while read line; do awk '{print $1}' $line | sort | uniq > $line.ids; done
rm *sam

# do same with rRNA hits
cd /scratch/amann3/rRNA_map
ls *sam | while read line; do awk '{print $1}' $line | sort | uniq > $line.ids; done
rm *sam

# make directory for filtered reads
cd /scratch/amann3
mkdir filtered

# concatenate filtered ids for each sample
cd /scratch/amann3/human_map
ls *ids | sed 's/.sam.ids//' | while read line; do cat $line.sam.ids /scratch/amann3/rRNA_map/$line.sam.ids | sort | uniq > /scratch/amann3/filtered/$line.filt.ids; done

echo 'Reads to be filtered concatenated'
```

Now merge together IDs, filter reads

```bash
# generate sample specific pbs scripts
cat /home/amann3/hiv_rnaseq/sample.ids | while read line; do sed "s/SAMPLE/$line/g" example.pbs > $line.pbs; done
```

Example PBS file

```bash
#PBS -N SAMPLE-remove_human_rRNA
#PBS -l select=1:ncpus=1:mem=750gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe
#PBS -q bigmem

# load modules
module add python/3.9.7

# now remove unwanted reads
cd /scratch/amann3/filtered
python /scratch/amann3/remove_seqs.py -f /scratch/amann3/merged_Qfiltered/SAMPLE.1.fastq.gz -i SAMPLE.filt.ids -o SAMPLE.1.fastq; done
python /scratch/amann3/remove_seqs.py -f /scratch/amann3/merged_Qfiltered/SAMPLE.2.fastq.gz -i SAMPLE.filt.ids -o SAMPLE.2.fastq; done

echo 'Filtering complete'
```

9. Denovo transcriptome assembly with RNASpades

```bash
#PBS -N SAMPLE-rnaSpades
#PBS -l select=1:ncpus=40:mem=750gb
#PBS -l walltime=138:00:00
#PBS -j oe
#PBS -m abe
#PBS -q bigmem

# change directory
cd /scratch/amann3/filtered/

# load module
module add spades/3.15.5

rnaspades.py -1 /scratch/amann3/filtered/SAMPLE.1.fastq -2 /scratch/amann3/filtered/SAMPLE.2.fastq -o /scratch/amann3/denovo_assembly/SAMPLE --checkpoints all -m 750 -t 40
```

