#!/bin/bash

#SBATCH --job-name denovo-genomeGenerate
#SBATCH --nodes 1
#SBATCH --tasks-per-node 10
#SBATCH --mem 750gb
#SBATCH --time 72:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# move to folder with data
cd /scratch/amann3/all_assemblies

# load module
module load star/2.7.5c

# run command
STAR --runMode genomeGenerate \
   --genomeFastaFiles /scratch/amann3/all_assemblies/all_assemblies.fna  \
   --runThreadN 10 \
   --limitGenomeGenerateRAM 228697571594 \
   --sjdbGTFfile /scratch/amann3/all_assemblies/all_assemblies.gtf \
   --genomeChrBinNbits 10 \
   --limitSjdbInsertNsj 3204994 \
   --genomeSAindexNbases 10\
   --sjdbGTFfeatureExon CDS
