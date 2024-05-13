# 1. Get the samples
```sh
conda activate 2024-hiv_rnaseq
cd ~/rna_dohmain/miseq_rna
cat samples | sed 's/_.*//' > ../rpoc/samples.txt
cd ../rpoc
ls ../../all_rpoc > rpoc_files
cat samples.txt | while read line; do grep $line ./rpoc_files; done > needed_rpoc
mkdir data
cd ../../all_rpoc
xargs -a ../rna_dohmain/rpoc/needed_rpoc cp -t ../rna_dohmain/rpoc/data
cd ~/rna_dohmain/rpoc/data
DM00428V2PQ84
```
# 2. Install packages
```R
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("dada2")
# install.packages("magrittr")
# install.packages("stringr")
# install.packages("data.table")
# install.packages("qualpalr")
# install.packages("seqinr")
```
# 3. Load packages
```R
library(dada2, verbose = FALSE)
library(stringr, verbose = FALSE)
library(data.table, verbose = FALSE)
library(qualpalr, verbose = FALSE)
library(ShortRead, verbose = FALSE)
library(Biostrings, verbose = FALSE)
library(seqinr, verbose = FALSE)
set.seed(12349)
```
# 4. Set up paths
```R
rawpath <- "/home/suzanne/rna_dohmain/rpoc/data"
fnFs <- sort(list.files(rawpath, pattern="_R1_001.fastq.gz", full.names=T))
fnRs <- sort(list.files(rawpath, pattern="_R2_001.fastq.gz", full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names, 50)
paste("Number of input samples: ", length(sample.names))
```
# 5. Make quality plot scores
```R
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
```
# 6. Preliminary read filter (remove any reads with Ns)
```R
fnFs.filtN <- file.path(rawpath, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(rawpath, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)
```
# 7. Cutadapt remove primer sequences and any adapter contamination
```R
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
```
# 8. Quality Filter and Trim Reads
```R
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, trimRight=25, maxN=c(0,0), maxEE=c(4,6), rm.phix=TRUE, matchIDs=TRUE, compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$precentage_retained <- retained$reads.out/retained$reads.in * 100
retained
```
# 9. Learn error rates
```R
errF <- learnErrors(filtFs, multithread=T, random=T)
errR <- learnErrors(filtRs, multithread=T, random=T)
err.f.plt <- plotErrors(errF, nominalQ=TRUE) 
err.f.plt
pdf(paste("img/", "error_plot.pdf", sep=""))
err.f.plt
dev.off()
system("~/.iterm2/imgcat ./img/error_plot.pdf")
```
# 10. Dereplicating
```R
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```
# 11. Sample infrence
```R
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, verbose = FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, verbose = FALSE)
```
# 12. Merge Sample reads
```R
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)
```
# 13. Generate sequence tables
```R
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```
# 14. Check Length distrubtion
```R
table(nchar(colnames(seqtab)))
```
# 15. Filter Sequences by length
```R
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
```
# 16. Remove chimeric sequences
```R
seqtab.nochim <- removeBimeraDenovo(seqtab.filt, method="pooled", multithread=T, verbose=T)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab.filt)
```
# 17. Processing summary
```R
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochimeras")
rownames(track) <- sample.names
track
# write to file
write.table(data.frame("row_names"=rownames(track),track),"read_retention.txt", row.names=FALSE, quote=F, sep="\t")
uniquesToFasta(seqtab.nochim, "rep_set.fa")
# fix ASV names 
system("awk '/^>/{print \">ASV\" ++i; next}{print}' < rep_set.fa > rep_set_fix.fa")
system("mv rep_set_fix.fa rep_set.fa")
# write sequence table to file, fix ASV names
my_otu_table <- t(as.data.frame(seqtab.nochim)) 
ASV.seq <- as.character(unclass(row.names(my_otu_table))) 
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') 
colnames(seqtab.nochim) <- ASV.num 
write.table(data.frame("row_names"=rownames(seqtab.nochim),seqtab.nochim),"sequence_table.merged.txt", row.names=FALSE, quote=F, sep="\t")
```
# 18. Assign taxonomy
```sh 
#taxonomic assignment
kraken2 --db ~/databases/kraken_rpoc_v0.2.0 \
	--threads 6 \
	--use-names \
	--output rep_set.kraken.out rep_set.fa \
	--unclassified-out rep_set.unclassified.out --confidence 0.01
```
# 19. Get full taxonomy
```sh
awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > taxids
sed "s/\t//g" ~/ref_db/ncbi_taxonomy/rankedlineage.dmp > rankedlineage.dmp
sort -t "|" -k 1b,1 rankedlineage.dmp > rankedlineage_sorted
cat rankedlineage_sorted | sed 's/|\{2,\}/|/g' > rankedlineage_clean
# add unclassified to taxonomy file
sed -i '1 i\0|unclassified|' rankedlineage_clean
```
# 20. Check if taxids exist in ranked lineage file
```sh
cat taxids | sed 's/$/|/' | while read line; do grep -c -m 1 ^$line rankedlineage_clean | sed "s/0/$line not found/"; done > missing_check
grep "not found" missing_check # should come back with nothing if all taxids found
cat taxids | sed 's/$/|/' | while read line; do grep -m 1 ^$line rankedlineage_clean || echo $line "no lineage" ; done > lineages #use the || pipe
grep "no lineage" lineages # should return empty
```
# 21. Prep taxonomy file
```sh
# remove taxids from file and reverse field order
sed 's/|/\t/' lineages | awk -F"\t" '{print $2}' | awk -F\| '{s=$NF;for(i=NF-1;i>=1;i--)s=s FS $i;print s}' | sed 's/^|//' | sed 's/ /_/g' | sed 's/|/;/g' > taxonomy
# merge asv ids and taxonomy
awk '{print $2}' rep_set.kraken.out > asvids
paste asvids taxonomy > taxonomy.txt
```
# 22. Finish taxonomy file
```sh
#filter out what was only assigned at phylum level
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt
```