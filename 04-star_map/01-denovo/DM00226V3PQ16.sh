#!/bin/bash

#SBATCH --job-name denovoMap-DM00226V3PQ16
#SBATCH --nodes 1
#SBATCH --tasks-per-node 40
#SBATCH --cpus-per-task 1
#SBATCH --mem 750gb
#SBATCH --time 72:00:00

#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# load module
module add biocontainers star/2.7.10b

# move into scratch
cd /scratch/amann3/filtered

# run command
STAR --runThreadN 40 \
   --genomeDir /scratch/amann3/all_assemblies/GenomeDir \
   --readFilesIn <(gunzip -c DM00226V3PQ16.1.fastq.gz) <(gunzip -c DM00226V3PQ16.2.fastq.gz) \
   --outFileNamePrefix /scratch/amann3/filtered/DM00226V3PQ16.denovo \
   --outSAMtype BAM Unsorted \
   --outReadsUnmapped /scratch/amann3/filtered/DM00226V3PQ16.denovo.unmapped.bam \
   --quantMode TranscriptomeSAM GeneCounts \
   --alignIntronMax 1 \
   --chimOutType SeparateSAMold
