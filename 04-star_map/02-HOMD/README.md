# HOMD mapping

```bash
#PBS -N homd_genome-generate
#PBS -l select=1:ncpus=25:mem=750gb
#PBS -l walltime=168:00:00
#PBS -j oe
#PBS -m abe
#PBS -q bigmem

# load module
module add biocontainers star/2.7.10b gffread/0.12.7

# move into scratch
cd /scratch/amann3/homd_db

# download data
wget https://www.homd.org/ftp/genomes/PROKKA/current/gff/ALL_genomes.gff
wget https://www.homd.org/ftp/genomes/PROKKA/current/fna/ALL_genomes.fna

# run command
gffread ALL_genomes.gff -T -o ALL_genomes.gtf

STAR --runMode genomeGenerate \
	--genomeFastaFiles ALL_genomes.fna  \
	--runThreadN 25 \
	--sjdbGTFfile ALL_genomes.gtf \
	--genomeChrBinNbits 10 \
	--limitGenomeGenerateRAM 1493556046773 \
	--sjdbGTFfeatureExon CDS
```

During this process, palmetto was updated to slurm, restart and rerun samples

```bash
#!/bin/bash

#SBATCH --job-name homdMap-SAMPLE
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
cd /scratch/amann3/filtered

# run command
STAR --runThreadN 40 \
	--genomeDir /scratch/amann3/homd_db/GenomeDir \
	--readFilesIn <(gunzip -c SAMPLE.1.fastq.gz) <(gunzip -c SAMPLE.2.fastq.gz) \
	--outFileNamePrefix /scratch/amann3/filtered/SAMPLE.homd \
	--outSAMtype BAM /scratch/amann3/filtered/Unsorted \
	--outReadsUnmapped /scratch/amann3/filtered/SAMPLE.homd.unmapped.bam \
	--quantMode TranscriptomeSAM GeneCounts \
	--alignIntronMax 1 \
	--chimOutType SeparateSAMold
```

Create new sample specific mapping jobs

```bash
cat /home/amann3/hiv_rnaseq/sample.ids | while read line; do sed "s/SAMPLE/$line/g" example.sh > $line.sh ; done
# submit
ls D*sh | while read line; do sbatch $line; done
```

Once the mapping is complete, download all output files to hillary and run feature counts

```bash
conda activate 2024-HIV_RNA
# conda install -c bioconda subread
cd /home/allie/domhain_RNAseq/03-star_map/02-HOMD_map
mkdir featurecounts
ls *Aligned.out.bam | sed 's/Aligned.out.bam//' | while read line; do featureCounts -f -p -C -B -M -a ../homd_db/ALL_genomes.gtf -o featurecounts/$line.out -T 60 $line\Aligned.out.bam -t CDS -g transcript_id; done
```

Combine output of feature counts into single file

```bash
cd featurecounts
paste *out | grep -v "^#" | awk '{printf "%s\t", $1}{for (i=7;i<=NF;i+=7) printf "%s\t", $i; printf "\n"}' > read_counts.txt
# clean up sample names
sed 's/Aligned.out.bam//g' read_counts.txt | sed 's/.homd//'g > temp
mv temp read_counts.txt
```










