#!/bin/bash

#SBATCH --job-name denovoMap-SAMPLE
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
   --genomeDir /scratch/amann3/denovo_db/GenomeDir \
   --readFilesIn <(gunzip -c SAMPLE.1.fastq.gz) <(gunzip -c SAMPLE.2.fastq.gz) \
   --outFileNamePrefix /scratch/amann3/filtered/SAMPLE.denovo \
   --outSAMtype BAM /scratch/amann3/filtered/Unsorted \
   --outReadsUnmapped /scratch/amann3/filtered/SAMPLE.denovo.unmapped.bam \
   --quantMode TranscriptomeSAM GeneCounts \
   --alignIntronMax 1 \
   --chimOutType SeparateSAMold
