#!/bin/bash

#SBATCH --job-name denovo-gtfgen
#SBATCH --nodes 1
#SBATCH --tasks-per-node 10
#SBATCH --mem 750gb
#SBATCH --time 72:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# move to folder with data
cd /scratch/amann3/all_assemblies

# load module
module load gffread/0.12.7 

# run
gffread -T all_assemblies.gff -o all_assemblies.gtf
echo "gffread complete"

