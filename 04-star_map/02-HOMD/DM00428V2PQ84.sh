#!/bin/bash

#SBATCH --job-name homdMap-DM00428V2PQ84
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
	--genomeDir /scratch/amann3/homd_db/GenomeDir \
	--readFilesIn <(gunzip -c DM00428V2PQ84.1.fastq.gz) <(gunzip -c DM00428V2PQ84.2.fastq.gz) \
	--outFileNamePrefix DM00428V2PQ84.homd \
	--outSAMtype BAM Unsorted \
	--outReadsUnmapped /scratch/amann3/filtered/DM00428V2PQ84.homd.unmapped.bam \
	--quantMode TranscriptomeSAM GeneCounts \
	--alignIntronMax 1 \
	--chimOutType SeparateSAMold
