# 1. 16S DADA2
```sh
mkdir usa_16S
cd usa_16S
mkdir data
cd data
sudo rsync /home/lauren/pugh_copy/16S_Plaque/raw_fastq/*fastq.gz suzanne@stella.clemson.edu:~/rna_dohmain/09-urease/00-scratch/usa_16S/data
scp Nigeria_USA_ITS1.tsv suzanne@stella.clemson.edu:~/rna_dohmain/09-urease/00-scratch/usa_16S/data
grep USA ~/usa_nigeria/Nigeria_USA_meta.txt > map.txt
head -n 1 ~/usa_nigeria/Nigeria_USA_meta.txt > header
cat header map.txt > temp
mv temp map.txt
#rename files
#R1
ls *_R1_*.fastq.gz | sort -n > num_order_R1.txt
grep UF Nigeria_USA_ITS1.tsv  | awk '{print $4}' > new_headers
paste new_headers <(sed 's/.*_S/_S/' num_order_R1.txt)  | sed 's/\t//'> new_name_R1.txt
paste -d"\t" num_order_R1.txt new_name_R1.txt > name_R1_map.txt
awk -F'/t' 'system("cp " $1 " " $2)' name_R1_map.txt
#R2
ls *_R2_*.fastq.gz | sort -n > num_order_R2.txt
grep UF Nigeria_USA_ITS1.tsv  | awk '{print $4}' > new_headers
paste new_headers <(sed 's/.*_S/_S/' num_order_R2.txt)  | sed 's/\t//' > new_name_R2.txt
paste -d"\t" num_order_R2.txt new_name_R2.txt > name_R2_map.txt
awk -F'/t' 'system("cp " $1 " " $2)' name_R2_map.txt
#move old files into different directory
mkdir old
cat num_order_R1.txt | while read line; do mv $line old; done
cat num_order_R2.txt | while read line; do mv $line old; done
```
# 2. DADA2 Pipeline
```R
#Load packages
library(dada2, verbose = FALSE)
library(stringr, verbose = FALSE)
library(data.table, verbose = FALSE)
library(ShortRead, verbose = FALSE)
library(Biostrings, verbose = FALSE)
library(seqinr, verbose = FALSE)
library(qualpalr, verbose = FALSE)
set.seed(12349)
#Set up paths
rawpath <- "/home/suzanne/rna_dohmain/09-urease/00-scratch/usa_16S/data"
fnFs <- sort(list.files(rawpath, pattern="_R1_001.fastq.gz", full.names=T))
fnRs <- sort(list.files(rawpath, pattern="_R2_001.fastq.gz", full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names, 50)
paste("Number of input samples: ", length(sample.names))
#Make quality plot scores
system("mkdir img") # ignore warning
fwdqual <- plotQualityProfile(fnFs[10:25])
revqual <- plotQualityProfile(fnRs[10:25])

pdf(paste("img/", "forward_quality_plot.pdf", sep=""))
fwdqual
dev.off()
pdf(paste("img/", "reverse_quality_plot.pdf", sep=""))
revqual
dev.off()
system("~/.iterm2/imgcat ./img/forward_quality_plot.pdf")
system("~/.iterm2/imgcat ./img/reverse_quality_plot.pdf")

#Filter out uncalled bases
fnFs.filtN <- file.path(rawpath, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(rawpath, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)
#Primer Removal
cutadapt <- as.character(system("which cutadapt", intern=T))
system("cutadapt --version")
path.cut <- file.path(rawpath, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc("GTGCCAGCMGCCGCGGTAA")
REV.RC <- dada2:::rc("GGACTACHVGGGTWTCTAAT")
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", "GTGCCAGCMGCCGCGGTAA", "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", "GGACTACHVGGGTWTCTAAT", "-A", FWD.RC)
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("--cores=0", R1.flags, R2.flags, "-n", 2,"-o", fnFs.cut[i], "-p", fnRs.cut[i], fnFs.filtN[i], fnRs.filtN[i]))
}
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2", full.names = TRUE))
#Quality Filter and Trim Reads
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, minLen = c(240,175),maxN=c(0,0), maxEE=c(2,2), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE,compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$precentage_retained <- retained$reads.out/retained$reads.in * 100
retained
# Learn error rates
errF <- learnErrors(filtFs, multithread=T, random=T)
errR <- learnErrors(filtRs, multithread=T, random=T)
err.f.plt <- plotErrors(errF, nominalQ=TRUE) 
pdf(paste("./img/", "error_plot.pdf", sep=""))
err.f.plt
dev.off()
system("~/.iterm2/imgcat ./img/error_plot.pdf")
# Dereplicating
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
# Sample infrence
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, verbose = FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, verbose = FALSE)
# Merge Sample reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)
# Generate sequence tables
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Check Length distrubtion 
table(nchar(colnames(seqtab)))
# get histogram of length distribution after filter
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab))))
len.plt <- plot(x=length.histogram[,1], y=length.histogram[,2])
pdf(paste("./img/", "length_hist.pdf", sep=""))
plot(x=length.histogram[,1], y=length.histogram[,2])
dev.off()
system("~/.iterm2/imgcat ./img/length_hist.pdf")
# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=T, verbose=T)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
# Processing summary
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochimeras")
rownames(track) <- sample.names
track
# write to file
write.table(data.frame("row_names"=rownames(track),track),"read_retention_16S.txt", row.names=FALSE, quote=F, sep="\t")
uniquesToFasta(seqtab.nochim, "rep_set.fa")
# fix ASV names 
system("awk '/^>/{print \">sASV\" ++i; next}{print}' < rep_set.fa > rep_set_fix.fa")
system("mv rep_set_fix.fa rep_set_16S.fa")
# write sequence table to file, fix ASV names
my_otu_table <- t(as.data.frame(seqtab.nochim)) 
ASV.seq <- as.character(unclass(row.names(my_otu_table))) 
ASV.num <- paste0("sASV", seq(ASV.seq), sep='') 
colnames(seqtab.nochim) <- ASV.num 
write.table(data.frame("row_names"=rownames(seqtab.nochim),seqtab.nochim),"sequence_table_16S.merged.txt", row.names=FALSE, quote=F, sep="\t")
save.image("16S_dada2.RData")
```
Make 16S database
```sh
mkdir database
cd database
grep ">" ../../../../homd_map/ALL_genomes.fna | sort | sed 's/|.*//' | uniq | sed 's/>//' > gff_files
grep "product=16S ribosomal RNA" ../../../../homd_map/ALL_genomes.gff | sed 's/\t/tab/g' | sed 's/.*ID=//g' | sed 's/;.*//' > 16S_ids
sed 's/_.*//' 16S_ids > seqIDs
#get seqID
wget https://www.homd.org/ftp/genomes/PROKKA/current/SEQID_info.txt
sed 's/ /\t/g' SEQID_info.txt | sed 's/\t.*https/\thttps/'| sed 's/\t.*GCA_/\tGCA_/' | sed 's/_[^_]*//2g' > seqIDs2GCA #format it so it is seqid and GCA
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
sed 's/ /_/g' assembly_summary_genbank.txt | awk '{print $1,$6}' | sed '1d' | sed '1d' | sed 's/ /\t/g' > gca2taxid #get just gca and taxid colum
sed 's/_.*//' 16S_ids > rpoc_seqs #format rpoc ref seqs with ID
for i in `cat rpoc_seqs`; do grep $i seqIDs2GCA | awk '{print $2}'; done > GCAs #get the GCA for refrence
parallel -a GCAs -j 7 -k "grep '{}' gca2taxid"> taxid
#check for missing ids
awk '{print $1}' taxid > tax_acc
awk '{print $2}' seqIDs2GCA | sed '1d'> GCA_rpoc 
cat GCA_rpoc tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids
parallel -a missed_ids -j 7 -k "grep '{}' gca2taxid"> taxid2
cat taxid taxid2 > all_taxids
awk '{print $1}' all_taxids > tax_acc
cat GCA_rpoc tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids
#find missing ID GCA
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank_historical.txt
for i in `cat missed_ids`; do grep -m 1 $i assembly_summary_genbank_historical.txt | awk '{print $1, $6}'; done > taxid3 #get the GCA for refrence
cat taxid taxid2 taxid3 > all_taxids
awk '{print $1}' all_taxids > tax_acc
cat GCA_rpoc tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids
#manually involving stuff
#BE CAREFUL HERE
sed 's/\..[1-2]*//' missed_ids| while read line; do grep -m 1 $line assembly_summary_genbank.txt | awk '{print $1, $6}'; done > taxid4 #get the GCA for refrence
sed  -i '1d' missed_ids
paste missed_ids taxid4 > temp
awk '{print $1, $3}' temp > taxid4
cat taxid taxid2 taxid3 taxid4 | sed 's/ /\t/' > all_taxids
awk '{print $1}' all_taxids > tax_acc
cat GCA_rpoc tax_acc > all_ids
grep -f <(sort all_ids | uniq -u) all_ids > missed_ids
#get just uniq GCA
sort all_taxids | uniq > GCA_2_taxid
sort rpoc_seqs | uniq | wc -l #sanity check
parallel -a rpoc_seqs -j 7 -k "grep '{}' seqIDs2GCA"> seqf2GCA_rpoc
python3 GCA2taxid.py

#fix the headers
paste 16S_ids taxids | sed 's/\t/|kraken:taxid|/' | sed 's/\t/ /' | sed 's/^/>/' > fixed_headers
#get sequences
cp ../../../02-operon-mapping/ALL_genomes.fna ./
cat 16S_ids | while read line; do grep -m 1 $line ../../../../homd_map/ALL_genomes.gff; done | awk -F "\t" '{print $1, $4, $5}'
grep "product=16S ribosomal RNA" ../../../../homd_map/ALL_genomes.gff > 16S.gff
parallel -a 16S_ids -j 7 -k "grep -wm 1 '{}' 16S.gff"| awk -F "\t" '{print $1, $4, $5}' | sed 's/ /\t/g' | sed 's/|/_/' > ids.txt
#subset the operon
# awk -v number="300" '$3!=0{$3+=number} 1' ids.txt > temp #adding 100 the end
awk -v number="1" '$3!=0{$3+=number} 1' ids.txt > temp #adding 1 the end
mv temp ids.txt
awk -v number="1" '$2!=0{$2+=number} 1' ids.txt > temp #adding 1 the end
mv temp ids.txt
seqtk subseq ALL_genomes.fna ids.txt > 16S_gens
grep -v ">" 16S_gens > seqs
paste fixed_headers seqs | sed 's/\t/\n/' > ref_16S.fa

mkdir kraken_homd
kraken2-build --download-taxonomy --db kraken_homd/
kraken2-build --add-to-library ref_16S.fa --db kraken_homd/
kraken2-build --build --max-db-size 8000000000 --db kraken_homd/
```
Assign taxonomy
```sh
#taxonomic assignment
kraken2 --db ./database/kraken_homd \
	--threads 6 \
	--use-names \
	--output rep_set.kraken.out rep_set_16S.fa \
	--unclassified-out rep_set.unclassified.out --confidence 0.01

awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > taxids
sed "s/\t//g" ../../../rpoc/database/rankedlineage.dmp > rankedlineage.dmp
sort -t "|" -k 1b,1 rankedlineage.dmp > rankedlineage_sorted
cat rankedlineage_sorted | sed 's/|\{2,\}/|/g' > rankedlineage_clean
# add unclassified to taxonomy file
sed -i '1 i\0|unclassified|' rankedlineage_clean
sed 's/|/\t/' rankedlineage_clean | sed 's/ /_/g' >rankedlineage_clean2
python3 lineages.py #output is lineage
awk -F"\t" '{print $2}' lineage | awk -F\| '{s=$NF;for(i=NF-1;i>=1;i--)s=s FS $i;print s}' | sed 's/^|//' | sed 's/ /_/g' | sed 's/|/;/g' > taxonomy
# merge assemebleies ids and taxonomy
awk '{print $2}' rep_set.kraken.out > asvids
paste asvids taxonomy > taxonomy.txt
#filter out what was only assigned at phylum level
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python ~/bin/fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt
```
# 2. rpoC DADA2
```sh
cat  ../../map.txt | awk '{print $1}' | sed '1d' | while read line; do ls ~/all_rpoc/$line*;done
```
DADA2
```R
#Load packages
library(dada2, verbose = FALSE)
library(stringr, verbose = FALSE)
library(data.table, verbose = FALSE)
library(ShortRead, verbose = FALSE)
library(Biostrings, verbose = FALSE)
library(seqinr, verbose = FALSE)
library(qualpalr, verbose = FALSE)
set.seed(12349)

#Set up paths
rawpath <- "~/rna_dohmain/09-urease/00-scratch/usa_rpoc/data"
fnFs <- sort(list.files(rawpath, pattern="_R1_001.fastq.gz", full.names=T))
fnRs <- sort(list.files(rawpath, pattern="_R2_001.fastq.gz", full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names, 50)
paste("Number of input samples: ", length(sample.names))
#Make quality plot scores
system("mkdir img") # ignore warning
fwdqual <- plotQualityProfile(fnFs[10:25])
revqual <- plotQualityProfile(fnRs[10:25])

pdf(paste("img/", "forward_quality_plot.pdf", sep=""))
fwdqual
dev.off()
pdf(paste("img/", "reverse_quality_plot.pdf", sep=""))
revqual
dev.off()
system("~/.iterm2/imgcat ./img/forward_quality_plot.pdf")
system("~/.iterm2/imgcat ./img/reverse_quality_plot.pdf")
#Filter out uncalled bases
fnFs.filtN <- file.path(rawpath, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(rawpath, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE, matchIDs = TRUE)
#Primer Removal
cutadapt <- as.character(system("which cutadapt", intern=T))
system("cutadapt --version")
path.cut <- file.path(rawpath, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc("MAYGARAARMGNATGYTNCARGA")
REV.RC <- dada2:::rc("GMCATYTGRTCNCCRTCRAA")
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", "MAYGARAARMGNATGYTNCARGA", "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", "GMCATYTGRTCNCCRTCRAA", "-A", FWD.RC)
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("--cores=0", R1.flags, R2.flags, "-n", 2,"-o", fnFs.cut[i], "-p", fnRs.cut[i], fnFs.filtN[i], fnRs.filtN[i]))
}
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2", full.names = TRUE))
# Quality Filter and Trim Reads
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, trimRight=25, maxN=c(0,0), maxEE=c(4,6), rm.phix=TRUE, matchIDs=TRUE, compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$precentage_retained <- retained$reads.out/retained$reads.in * 100
retained
# Learn error rates
errF <- learnErrors(filtFs, multithread=T, random=T)
errR <- learnErrors(filtRs, multithread=T, random=T)
err.f.plt <- plotErrors(errF, nominalQ=TRUE) 
err.f.plt
pdf(paste("img/", "error_plot.pdf", sep=""))
err.f.plt
dev.off()
system("~/.iterm2/imgcat ./img/error_plot.pdf")
# Derepelicating
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
# Sample infrence
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, verbose = FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, verbose = FALSE)
# Merge Reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)
# Generate sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Check length distro
table(nchar(colnames(seqtab)))
# Filter sequence by length
seqlens <- nchar(getSequences(seqtab))
seqtab.filt <- seqtab[,seqlens >= 450]
dim(seqtab.filt)
# get histogram of length distribution after filter
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab.filt))))
len.plt <- plot(x=length.histogram[,1], y=length.histogram[,2])
pdf(paste("img/", "length_hist.pdf", sep=""))
plot(x=length.histogram[,1], y=length.histogram[,2])
dev.off()
system("~/.iterm2/imgcat ./img/length_hist.pdf")
# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab.filt, method="pooled", multithread=T, verbose=T)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab.filt)
# Processing summary
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochimeras")
rownames(track) <- sample.names
track
# write to file
write.table(data.frame("row_names"=rownames(track),track),"read_retention.txt", row.names=FALSE, quote=F, sep="\t")
uniquesToFasta(seqtab.nochim, "rep_set.fa")
# fix ASV names 
system("awk '/^>/{print \">rASV\" ++i; next}{print}' < rep_set.fa > rep_set_fix.fa")
system("mv rep_set_fix.fa rep_set.fa")
# write sequence table to file, fix ASV names
my_otu_table <- t(as.data.frame(seqtab.nochim)) 
ASV.seq <- as.character(unclass(row.names(my_otu_table))) 
ASV.num <- paste0("rASV", seq(ASV.seq), sep='') 
colnames(seqtab.nochim) <- ASV.num 
write.table(data.frame("row_names"=rownames(seqtab.nochim),seqtab.nochim),"sequence_table.merged.txt", row.names=FALSE, quote=F, sep="\t")
save.image("rpoC_dada2.RData")
```
Assign taxonomy
```sh
kraken2 --db ~/rna_dohmain/rpoc/database/kraken_homd \
	--threads 6 \
	--use-names \
	--output rep_set.kraken.out rep_set.fa \
	--unclassified-out rep_set.unclassified.out --confidence 0.01
awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > taxids
sed "s/\t//g" ../../../rpoc/database/rankedlineage.dmp > rankedlineage.dmp
sort -t "|" -k 1b,1 rankedlineage.dmp > rankedlineage_sorted
cat rankedlineage_sorted | sed 's/|\{2,\}/|/g' > rankedlineage_clean
# add unclassified to taxonomy file
sed -i '1 i\0|unclassified|' rankedlineage_clean
sed 's/|/\t/' rankedlineage_clean | sed 's/ /_/g' >rankedlineage_clean2
python3 ../usa_16S/lineages.py #output is lineage
awk -F"\t" '{print $2}' lineage | awk -F\| '{s=$NF;for(i=NF-1;i>=1;i--)s=s FS $i;print s}' | sed 's/^|//' | sed 's/ /_/g' | sed 's/|/;/g' > taxonomy
# merge assemebleies ids and taxonomy
awk '{print $2}' rep_set.kraken.out > asvids
paste asvids taxonomy > taxonomy.txt
#filter out what was only assigned at phylum level
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python ~/bin/fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt
```
# 3. Compare distrubtion of Actionmyces genus
16S
```R
library(philr, warn.conflicts = F, quietly = T)
library(RColorBrewer, warn.conflicts = F, quietly = T)
library(UpSetR, warn.conflicts = F, quietly = T)
library(ggfortify, warn.conflicts = F, quietly = T)
library(randomForest, warn.conflicts = F, quietly = T)
# library(rfUtilities, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(gridExtra, warn.conflicts = F, quietly = T)
library(microbiome, warn.conflicts = F, quietly = T)
library(phylofactor, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(pairwiseAdonis, warn.conflicts = F, quietly = T)
library(ape, warn.conflicts = F, quietly = T)
library(metagMisc, warn.conflicts = F, quietly = T)
library(ranacapa, warn.conflicts = F, quietly = T)
library(MASS, warn.conflicts = F, quietly = T)
library(ggdendro, warn.conflicts = F, quietly = T)
#set working directory
setwd("~/rna_dohmain/09-urease/00-scratch")
#convert to relative abundance
seqtab <- t(read.table("~/rna_dohmain/09-urease/00-scratch/usa_16S/sequence_table_16S.merged.txt", header=T, row.names=1))
tax <- read.table("~/rna_dohmain/09-urease/00-scratch/usa_16S/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("~/rna_dohmain/09-urease/00-scratch/map.txt", sep="\t", header=T, row.names=1)
notinmeta <- setdiff(colnames(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), colnames(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw

ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
ps.dat

rel <- microbiome::transform(ps.dat, "compositional")
actino <- subset_taxa(rel, V7=="Actinomyces")
glom <- tax_glom(actino, taxrank=rank_names(actino)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
# plot
pdf("actino_16S.pdf", width =15, heigh =10)
ggplot(data)+
  geom_bar(aes(x=Sample, y=Abundance,fill=V8),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  	  axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
  	  legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./actino_16S.pdf")
test <- data %>% group_by(Sample) %>% summarise(average_baseMean = across(c(Abundance)))
glom <- tax_glom(actino, taxrank=rank_names(actino)[6])
data <- psmelt(glom) # create dataframe from phyloseq object
mean(data$Abundance)
````
rpoC
```R
library(philr, warn.conflicts = F, quietly = T)
library(RColorBrewer, warn.conflicts = F, quietly = T)
library(UpSetR, warn.conflicts = F, quietly = T)
library(ggfortify, warn.conflicts = F, quietly = T)
library(randomForest, warn.conflicts = F, quietly = T)
# library(rfUtilities, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(gridExtra, warn.conflicts = F, quietly = T)
library(microbiome, warn.conflicts = F, quietly = T)
library(phylofactor, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(pairwiseAdonis, warn.conflicts = F, quietly = T)
library(ape, warn.conflicts = F, quietly = T)
library(metagMisc, warn.conflicts = F, quietly = T)
library(ranacapa, warn.conflicts = F, quietly = T)
library(MASS, warn.conflicts = F, quietly = T)
library(ggdendro, warn.conflicts = F, quietly = T)
#set working directory
setwd("~/rna_dohmain/09-urease/00-scratch")
#convert to relative abundance
seqtab <- t(read.table("~/rna_dohmain/09-urease/00-scratch/usa_rpoc/sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("~/rna_dohmain/09-urease/00-scratch/usa_rpoc/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
map <- read.table("~/rna_dohmain/09-urease/00-scratch/map.txt", sep="\t", header=T, row.names=1)
notinmeta <- setdiff(colnames(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), colnames(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw

ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
ps.dat

rel <- microbiome::transform(ps.dat, "compositional")
actino <- subset_taxa(rel, V7=="Actinomyces")
glom <- tax_glom(actino, taxrank=rank_names(actino)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
data
data$Sample<- factor(data$Sample,levels=unique(data$Sample))
# plot
pdf("actino_rpoc.pdf", width =15, heigh =10)
ggplot(data)+
  geom_bar(aes(x=Sample, y=Abundance,fill=V8),stat="identity", position="stack")+
  facet_grid(~ tooth_health, switch = "x", scales = "free_x") +
  theme_minimal()+
  ggplot2::theme(axis.text.x = element_text(color = "black", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"),
  	  axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
  	  legend.text = element_text(size = 7),legend.title = element_text(size = 10))
dev.off()
system("~/.iterm2/imgcat ./actino_rpoc.pdf")
test <- data %>% group_by(Sample) %>% summarise(average_baseMean = across(c(Abundance)))
glom <- tax_glom(actino, taxrank=rank_names(actino)[6])
data <- psmelt(glom) # create dataframe from phyloseq object
mean(data$Abundance)
````