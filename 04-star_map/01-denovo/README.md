### Star mapping to denovo assembled data

First need to merge all of our denovo assemblies into one, dereplicate, and run prokka

```bash
#PBS -N prokka-denovo
#PBS -l select=1:ncpus=16:mem=100gb
#PBS -l walltime=72:00:00
#PBS -j oe

# move to folder with data
cd /fastscratch/amann3/denovo_assembly

# load modules
module add prokka/1.14.5 vsearch/2.21.1 seqtk/1.3-r106

# merge and dereplicate
cat */transcripts.fasta > all_assemblies.fasta
vsearch --derep_fulllength all_assemblies.fasta \
   --output all_assemblies.uniq.fasta
echo "all transcripts"
grep ">" all_assemblies.fasta
echo "unique transcripts"
grep ">" all_assemblies.uniq.fasta

# rename contigs so they are compliant with prokka (less than or equal to 37 characters long)
seqtk rename all_assemblies.uniq.fasta denovo_ > all_assemblies.shortid.fasta

# run prokka
prokka --cpus 16 --force --prefix all_assemblies  all_assemblies.shortid.fasta
echo "prokka complete"
```

Running out of walltime on prokka for some reason, try running on hillary instead to debug

```bash
cd /home/allie/domhain_RNAseq/03-star_map
# conda create -n 2024-HIV_RNA
conda activate 2024-HIV_RNA
# update conda
conda install conda=24.3.0
# conda install -y -c biobuilds perl=5.22
# conda install -y -c conda-forge parallel
conda install -y -c bioconda prodigal 
conda install -y -c bioconda blast=2.2
conda install -y -c bioconda tbl2asn
conda install -y -c bioconda prokka

prokka --cpus 60 \
   --force \
   --prefix all_assemblies  \
   --norrna \
   --notrna \
   all_assemblies.shortid.fasta \
   1> prokka.out \
   2> prokka.err
```

This is getting stuck on the tbl2asn step -- since I probably won't need a genbank file, run the rest of the scripts to get mapping.

Generate reference database from these results for star (now using slurm)

```bash
#!/bin/bash

#SBATCH --job-name denovo-genomeGenerate
#SBATCH --nodes 1
#SBATCH --tasks-per-node 10
#SBATCH --mem 750gb
#SBATCH --time 72:00:00r
#SBATCH --constraint interconnect_fdr

# move to folder with data
cd /scratch/amann3/all_assemblies

# load module
module add gffread/0.12.7 star/2.7.5c

# run
gffread -T all_assemblies.gff -o all_assemblies.gtf
echo "gffread complete"

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
```

Map filtered reads to denovo assembled contig database

```bash
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
```

Create new sample specific mapping jobs

```bash
cat /home/amann3/hiv_rnaseq/sample.ids | while read line; do sed "s/SAMPLE/$line/g" example.sh > $line.sh ; done
# submit
ls D*sh | while read line; do sbatch $line; done
```

Once mapping is complete, download all output files to hillary and run feature counts



```bash
conda activate 2024-HIV_RNA
# conda install -c bioconda subread
cd /home/allie/domhain_RNAseq/03-star_map/01-denovo
mkdir featurecounts
ls *Aligned.out.bam | sed 's/Aligned.out.bam//' | while read line; do featureCounts -f -p -C -B -M -a ../homd_db/ALL_genomes.gtf -o featurecounts/$line.out -T 60 $line\Aligned.out.bam -t CDS -g transcript_id; done
```

Combine output of feature counts into single file






paste *out | grep -v "^#" | awk '{printf "%s\t", $1}{for (i=7;i<=NF;i+=7) printf "%s\t", $i; printf "\n"}' > read_counts.txt
# sanity check that the genes you are looking for are in the file
grep "arcA" ../ALL_genomes.gtf | awk -F";" '{print $1}' | awk -F"\t" '{print $9}' | sed 's/transcript_id "//' | sed 's/"//' | sort | uniq > test_arcA.ids
cat test_arcA.ids | while read line; do grep $line read_counts.txt ; done
# clean up sample names
sed 's/Aligned.out.bam//g' read_counts.txt | sed 's/_S..\t/\t/'g | sed 's/_S.\t/\t/g' > temp
mv temp read_counts.txt










