# rpoC amplicon processing

### 1. Load environment

```bash
cd ~/domhain_RNAseq/11-rpoc_processing
conda activate 2024-HIV_RNASeq
conda install bioconda::cutadapt=
```

### 2. Install R packages

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")
install.packages("magrittr")
install.packages("stringr")
install.packages("data.table")
install.packages("broom")
install.packages("qualpalr")
install.packages("seqinr")
```

### 3. Load required R libraries

```R
library(dada2, warn.conflicts = F, quietly = T)
library(stringr, warn.conflicts = F, quietly = T)
library(data.table, warn.conflicts = F, quietly = T)
library(qualpalr, warn.conflicts = F, quietly = T)
library(ShortRead, warn.conflicts = F, quietly = T)
library(Biostrings, warn.conflicts = F, quietly = T)
library(seqinr, warn.conflicts = F, quietly = T)
sessionInfo()
```     

### 4. File path setup 

```R
rawpath <- "raw"
wdpath <- "~/domhain_RNAseq/11-rpoc_processing" # change to where git repository was cloned
fnFs <- sort(list.files(rawpath, pattern="_R1_001.fastq.gz", full.names=T))
fnRs <- sort(list.files(rawpath, pattern="_R2_001.fastq.gz", full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names, 10)
paste("Number of input samples: ", length(sample.names))
```

### 5. Plot quality scores

```R
fwdqual <- plotQualityProfile(fnFs[10:25])
revqual <- plotQualityProfile(fnRs[10:25])

pdf("forward_quality_plot.pdf")
fwdqual
dev.off()
pdf("reverse_quality_plot.pdf")
revqual
dev.off()
system("/home/allie/.iterm2/imgcat forward_quality_plot.pdf") # this is so I can view within the terminal, will need imgcat to run
system("/home/allie/.iterm2/imgcat reverse_quality_plot.pdf")
```

### 6. Preliminary filter (removes sequences with uncalled bases)

```R
fnFs.filtN <- file.path(rawpath, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(rawpath, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)
```

### 7. Primer removal with cutadapt

```R
cutadapt <- as.character("/home/allie/miniforge3/bin/cutadapt")
cutadapt
system(paste(cutadapt, "--version", sep=" "))
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

### 8. Final filter and trim reads

```R
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, trimRight=25, maxN=c(0,0), maxEE=c(4,6), rm.phix=TRUE, matchIDs=TRUE, compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
retained # what proportion of reads retained after filtering
```

### 9. Learn and plot error rates

```R
set.seed(12349)
errF <- learnErrors(filtFs, multithread=T, random=T)
errR <- learnErrors(filtRs, multithread=T, random=T)
err.f.plt <- plotErrors(errF, nominalQ=TRUE) 
pdf("error_plot.pdf")
err.f.plt
dev.off()
system("/home/allie/.iterm2/imgcat error_plot.pdf")
```
  
### 10. Dereplication

```R
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# reassign sample names
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

### 11. Sample inference

```R
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```     

### 12. Merge paired-end reads

```R
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=T)
```

### 13. Construct sequence table

```R
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

### 14. Length filter

```R
# first what does our length distribution look like?
table(nchar(colnames(seqtab)))  
seqlens <- nchar(getSequences(seqtab))
seqtab.filt <- seqtab[,seqlens >= 450] # only removing those below 450bp, from HOMD range should be around 477-493
dim(seqtab.filt)
# sequence length post filter
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab.filt))))
plot(x=length.histogram[,1], y=length.histogram[,2])
pdf("length_hist.pdf")
plot(x=length.histogram[,1], y=length.histogram[,2])
dev.off()
system("/home/allie/.iterm2/imgcat length_hist.pdf")
```

### 15. Remove chimeric sequences

```R
seqtab.nochim <- removeBimeraDenovo(seqtab.filt, method="pooled", multithread=T, verbose=T)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab.filt)
```

### 16. Processing summary

```R
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochimeras")
rownames(track) <- sample.names
track
```

### 17. Save output

```R
write.table(data.frame("row_names"=rownames(track),track),"read_retention.txt", row.names=FALSE, quote=F, sep="\t")
uniquesToFasta(seqtab.nochim, "rep_set.fa")
system("awk '/^>/{print \">ASV\" ++i; next}{print}' < rep_set.fa > rep_set_fix.fa")
system("mv rep_set_fix.fa rep_set.fa")
```

### 18. Clean up ASV names and write sequence table to file

```R
my_otu_table <- t(as.data.frame(seqtab.nochim)) 
ASV.seq <- as.character(unclass(row.names(my_otu_table))) 
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') 
colnames(seqtab.nochim) <- ASV.num 
write.table(data.frame("row_names"=rownames(seqtab.nochim),seqtab.nochim),"sequence_table.merged.txt", row.names=FALSE, quote=F, sep="\t")
```



     








