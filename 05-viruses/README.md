Quickly identify potential viral contigs from our denovo assemblies for Eduardo

```bash
cd /home/allie/domhain_RNAseq/02-denovo_assembly
cat */transcripts.fasta > all_transcripts.fasta
# dereplicate
vsearch --threads 25 --sizeout --derep_fulllength all_transcripts.fasta --output all_transcripts.uniq.fa
# download virus refseq database
mkdir ~/domhain_RNAseq/04-viruses && cd ~/domhain_RNAseq/04-viruses
mkdir k2_refviral && cd k2_refviral
wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz
tar xzf k2_viral_20240112.tar.gz
rm k2_viral_20240112.tar.gz
cd ..
# run kraken2
kraken2 --db k2_refviral \
	--threads 25 \
	--use-names \
	--classified-out refviral.classified.out \
	--unclassifed-out refviral.unclassified.out \
	--output refviral.kraken.out \
	../02-denovo_assembly/all_transcripts.uniq.fa \
	--confidence 0.01
```
