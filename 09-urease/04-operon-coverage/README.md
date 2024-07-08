# 1. Coverage plot for all upregulated genes
Grab the sequences 
```sh
#upload to palmetto
rsync -a ~/rna_dohmain/10-urease/02-operon-mapping/ure_operon.shortid.fasta scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/operon-coverage/
```
# 2. Map to operon on Palmetto
Build database
```sh
module add python/3.4.2 bowtie2/2.3.5.1
bowtie2-build ure_operon.shortid.fasta ureABC_operon.db

```
Run mapping
```sh
#!/bin/bash

#SBATCH --job-name operoncoverage-SAMPLE
#SBATCH --nodes 1
#SBATCH --tasks-per-node 10
#SBATCH --cpus-per-task 1
#SBATCH --mem 500gb
#SBATCH --time 72:00:00
#SBATCH --constraint interconnect_fdr
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# load module
module add samtools/1.12  python/3.4.2 bowtie2/2.3.5.1
# move into scratch
cd /scratch/scrull/hiv_rnaseq/operon-coverage

#map to ureABC
bowtie2 -x ./ureABC_operon.db -1 ../filtered/SAMPLE.1.fastq.gz -2 ../filtered/SAMPLE.2.fastq.gz --threads 10 --end-to-end  --qc-filter --no-unal -t -S ./SAMPLE.ureABC.sam 2>./SAMPLE.ureABC.out 1>./SAMPLE.ureABC.err
# convert to bam and sort by coordinate
samtools sort SAMPLE.ureABC.sam -o SAMPLE.ureABC.sort.bam
```
Create new sample specific mapping
```sh
ls /scratch/scrull/hiv_rnaseq/filtered | sed 's/\..*//' | sort | uniq | while read line; do sed "s/SAMPLE/$line/g" example.sh > $line.sh ; done
# submit
ls D*sh | while read line; do sbatch $line; done
```
Merge bam files by health status
```sh
# H 
grep -w H map.txt | awk '{print $1}' | while read line; do ls *bam| grep $line ; done > all_H.txt
samtools merge all_H.sort.bam -b all_H.txt
# D
grep -w D map.txt | awk '{print $1}' | while read line; do ls *bam| grep $line ; done > all_D.txt
samtools merge all_D.sort.bam -b all_D.txt
ls all*bam | while read line; do samtools depth $line > $line.depth; done

#HUU
grep -w H map.txt | grep -w HUU | awk '{print $1}' | while read line; do ls *bam| grep $line ; done > HUU_H.txt
samtools merge HUU_H.sort.bam -b HUU_H.txt
# D
grep -w D map.txt | grep -w HUU | awk '{print $1}' | while read line; do ls *bam| grep $line ; done > HUU_D.txt
samtools merge HUU_D.sort.bam -b HUU_D.txt
ls HUU*bam | while read line; do samtools depth $line > $line.depth; done
#HEU
grep -w H map.txt | grep -w HEU | awk '{print $1}' | while read line; do ls *bam| grep $line ; done > HEU_H.txt
samtools merge HEU_H.sort.bam -b HEU_H.txt
# D
grep -w D map.txt | grep -w HEU | awk '{print $1}' | while read line; do ls *bam| grep $line ; done > HEU_D.txt
samtools merge HEU_D.sort.bam -b HEU_D.txt
ls HEU*bam | while read line; do samtools depth $line > $line.depth; done
#HI
grep -w H map.txt | grep -w HI | awk '{print $1}' | while read line; do ls *bam| grep $line ; done > HI_H.txt
samtools merge HI_H.sort.bam -b HI_H.txt
# D
grep -w D map.txt | grep -w HI | awk '{print $1}' | while read line; do ls *bam| grep $line ; done > HI_D.txt
samtools merge HI_D.sort.bam -b HI_D.txt
ls HI*bam | while read line; do samtools depth $line > $line.depth; done
```
Transfer files to stella
```sh
scp 'scrull@slogin.palmetto.clemson.edu:/scratch/scrull/hiv_rnaseq/operon-coverage/*depth' ~/rna_dohmain/10-urease/04-operon-coverage
```
# 3. Make R plots for all depth
```R
library(ggplot2)
library(tidyverse)
setwd("~/rna_dohmain/10-urease/04-operon-coverage")
pd <- read.table("all_D.sort.bam.depth", header=F, sep="\t")
pf <- read.table("all_H.sort.bam.depth", header=F, sep="\t")
all.equal(pd,pf)

#get DESeq results
all_res <- read.table("../03-diff-abundance/deseq_results_operon-HvD.txt", header=T, sep="\t")
#get mapping file
map <- read.table("../02-operon-mapping/new.cords2", header=T, sep="\t")
annots <- read.table("../02-operon-mapping/operon_tre.annots", header=T, sep="\t")
annots$SEQ <- gsub(x = annots$SEQ, pattern = "-", replacement = ":")
annots$SEQ <- gsub(x = annots$SEQ, pattern = "_R_", replacement = "")
map <- inner_join(map, annots, by = join_by(full.1 == SEQ), keep=NULL)
#combine files
# get ids
all_res$seqs <- gsub(x = rownames(all_res), pattern = "\\..*", replacement = "")
#merge
all_res$ids <- rownames(all_res)
res_A <- inner_join(all_res, map, by = join_by(ids == idA), keep=TRUE)
res_B <- inner_join(all_res, map, by = join_by(ids == idB), keep=TRUE)
res_C <- inner_join(all_res, map, by = join_by(ids == idC), keep=TRUE)
resdf <- rbind(res_A, res_B, res_C)
# get list of loci that are significant and have a log FC of 2+
# set LFC and P value filters
lfc = 2
pval = 0.05
sig_df <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 
# pd <- rename(pd, seq = V1, pos = V2, depth =V3, A1=V4, A2 =V5, B1=V6, B2=V7, C1= V8, C2= V9)
pf <- inner_join(pf, map, by = join_by(V1 == full.1), keep=NULL)
seqs <- unique(sort(sig_df$full.1))
for(i in seqs) {
sub_H <- subset(pf, V1 == i )
sub_D <- subset(pd, V1 == i )

# print(all.equal(sub_H,sub_D))
df <- resdf %>% filter(full.1 == i)
x <- unique(sort(df$taxa))

y <- gsub(x=i, pattern = "_.*", replacement="")
p <- ggplot() +
  geom_line(data = sub_H, aes(x=V2, y=log10(V3)), color = "#44CCD0")+
  geom_line(data = sub_D, aes(x=V2, y=log10(V3)), color="#FA918B")+
  scale_x_continuous(name="Base Position") +
  scale_y_continuous(name="log10 coverage")+
  ggtitle(x)+
  geom_vline(xintercept = sub_H$newA1, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newA2, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newB1, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newB2, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newC1, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newC2, color = "black", size=.2)+
  theme_minimal()
pdf(paste0(y, "_cov",".pdf"))
print(p)
dev.off()
}

system("ls | grep pdf | grep -v H | while read line; do ~/.iterm2/imgcat $line;done")
ureABC_name <- names(which(table(sig_df$seqs) >= 3))
```
# 4. Make HUU Coverage plots
```R
library(ggplot2)
library(tidyverse)
pd <- read.table("HUU_D.sort.bam.depth", header=F, sep="\t")
pf <- read.table("HUU_H.sort.bam.depth", header=F, sep="\t")
all.equal(pd,pf)

#get DESeq results
all_res <- read.table("../03-diff-abundance/deseq_results_operon-HvD-HUU.txt", header=T, sep="\t")
#get mapping file
map <- read.table("../02-operon-mapping/new.cords2", header=T, sep="\t")
annots <- read.table("../02-operon-mapping/operon_tre.annots", header=T, sep="\t")
annots$SEQ <- gsub(x = annots$SEQ, pattern = "-", replacement = ":")
annots$SEQ <- gsub(x = annots$SEQ, pattern = "_R_", replacement = "")
map <- inner_join(map, annots, by = join_by(full.1 == SEQ), keep=NULL)
#combine files
# get ids
all_res$seqs <- gsub(x = rownames(all_res), pattern = "\\..*", replacement = "")
#merge
all_res$ids <- rownames(all_res)
res_A <- inner_join(all_res, map, by = join_by(ids == idA), keep=TRUE)
res_B <- inner_join(all_res, map, by = join_by(ids == idB), keep=TRUE)
res_C <- inner_join(all_res, map, by = join_by(ids == idC), keep=TRUE)
resdf <- rbind(res_A, res_B, res_C)
# get list of loci that are significant and have a log FC of 2+
# set LFC and P value filters
lfc = 2
pval = 0.05
sig_df <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# pd <- rename(pd, seq = V1, pos = V2, depth =V3, A1=V4, A2 =V5, B1=V6, B2=V7, C1= V8, C2= V9)
pf <- inner_join(pf, map, by = join_by(V1 == full.1), keep=NULL)
seqs <- unique(sort(sig_df$full.1))
for(i in seqs) {
sub_H <- subset(pf, V1 == i )
sub_D <- subset(pd, V1 == i )

# print(all.equal(sub_H,sub_D))
df <- resdf %>% filter(full.1 == i)
x <- unique(sort(df$taxa))

y <- gsub(x=i, pattern = "_.*", replacement="")
p <- ggplot() +
  geom_line(data = sub_H, aes(x=V2, y=log10(V3)), color = "#44CCD0")+
  geom_line(data = sub_D, aes(x=V2, y=log10(V3)), color="#FA918B")+
  scale_x_continuous(name="Base Position") +
  scale_y_continuous(name="log10 coverage")+
  ggtitle(x)+
  geom_vline(xintercept = sub_H$newA1, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newA2, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newB1, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newB2, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newC1, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newC2, color = "black", size=.2)+
  theme_minimal()
pdf(paste0(y, "_HUU_cov",".pdf"))
print(p)
dev.off()
}

system("ls | grep _HUU_cov | while read line; do ~/.iterm2/imgcat $line;done")
ureABC_HUU_name <- names(which(table(sig_df$seqs) >= 3))

```
# 5. Make HEU Coverage plots
```R
library(ggplot2)
library(tidyverse)
pd <- read.table("HEU_D.sort.bam.depth", header=F, sep="\t")
pf <- read.table("HEU_H.sort.bam.depth", header=F, sep="\t")
all.equal(pd,pf)

#get DESeq results
all_res <- read.table("../03-diff-abundance/deseq_results_operon-HvD-HEU.txt", header=T, sep="\t")
#get mapping file
map <- read.table("../02-operon-mapping/new.cords2", header=T, sep="\t")
annots <- read.table("../02-operon-mapping/operon_tre.annots", header=T, sep="\t")
annots$SEQ <- gsub(x = annots$SEQ, pattern = "-", replacement = ":")
annots$SEQ <- gsub(x = annots$SEQ, pattern = "_R_", replacement = "")
map <- inner_join(map, annots, by = join_by(full.1 == SEQ), keep=NULL)
#combine files
# get ids
all_res$seqs <- gsub(x = rownames(all_res), pattern = "\\..*", replacement = "")
#merge
all_res$ids <- rownames(all_res)
res_A <- inner_join(all_res, map, by = join_by(ids == idA), keep=TRUE)
res_B <- inner_join(all_res, map, by = join_by(ids == idB), keep=TRUE)
res_C <- inner_join(all_res, map, by = join_by(ids == idC), keep=TRUE)
resdf <- rbind(res_A, res_B, res_C)
# get list of loci that are significant and have a log FC of 2+
# set LFC and P value filters
lfc = 2
pval = 0.05
sig_df <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# pd <- rename(pd, seq = V1, pos = V2, depth =V3, A1=V4, A2 =V5, B1=V6, B2=V7, C1= V8, C2= V9)
pf <- dplyr::inner_join(pf, map, by = join_by(V1 == full.1), keep=NULL)
seqs <- unique(sort(sig_df$full.1))
for(i in seqs) {
sub_H <- subset(pf, V1 == i )
sub_D <- subset(pd, V1 == i )

# print(all.equal(sub_H,sub_D))
df <- resdf %>% filter(full.1 == i)
x <- unique(sort(df$taxa))

y <- gsub(x=i, pattern = "_.*", replacement="")
p <- ggplot() +
  geom_line(data = sub_H, aes(x=V2, y=log10(V3)), color = "#44CCD0")+
  geom_line(data = sub_D, aes(x=V2, y=log10(V3)), color="#FA918B")+
  scale_x_continuous(name="Base Position") +
  scale_y_continuous(name="log10 coverage")+
  ggtitle(x)+
  geom_vline(xintercept = sub_H$newA1, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newA2, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newB1, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newB2, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newC1, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newC2, color = "black", size=.2)+
  theme_minimal()
pdf(paste0(y, "_HEU_cov",".pdf"))
print(p)
dev.off()
}

system("ls | grep _HEU_cov | while read line; do ~/.iterm2/imgcat $line;done")
ureABC_HEU_name <- names(which(table(sig_df$seqs) >= 3))

```
# 6. Make HI Coverage plots
```R
library(ggplot2)
library(tidyverse)
pd <- read.table("HI_D.sort.bam.depth", header=F, sep="\t")
pf <- read.table("HI_H.sort.bam.depth", header=F, sep="\t")
all.equal(pd,pf)

#get DESeq results
all_res <- read.table("../03-diff-abundance/deseq_results_operon-HvD-HI.txt", header=T, sep="\t")
#get mapping file
map <- read.table("../02-operon-mapping/new.cords2", header=T, sep="\t")
annots <- read.table("../02-operon-mapping/operon_tre.annots", header=T, sep="\t")
annots$SEQ <- gsub(x = annots$SEQ, pattern = "-", replacement = ":")
annots$SEQ <- gsub(x = annots$SEQ, pattern = "_R_", replacement = "")
map <- inner_join(map, annots, by = join_by(full.1 == SEQ), keep=NULL)
#combine files
# get ids
all_res$seqs <- gsub(x = rownames(all_res), pattern = "\\..*", replacement = "")
#merge
all_res$ids <- rownames(all_res)
res_A <- inner_join(all_res, map, by = join_by(ids == idA), keep=TRUE)
res_B <- inner_join(all_res, map, by = join_by(ids == idB), keep=TRUE)
res_C <- inner_join(all_res, map, by = join_by(ids == idC), keep=TRUE)
resdf <- rbind(res_A, res_B, res_C)
# get list of loci that are significant and have a log FC of 2+
# set LFC and P value filters
lfc = 2
pval = 0.05
sig_df <- resdf %>% filter(padj <= pval) %>% filter(log2FoldChange >= lfc | log2FoldChange <= -lfc) 

# pd <- rename(pd, seq = V1, pos = V2, depth =V3, A1=V4, A2 =V5, B1=V6, B2=V7, C1= V8, C2= V9)
pf <- inner_join(pf, map, by = join_by(V1 == full.1), keep=NULL)
seqs <- unique(sort(sig_df$full.1))
for(i in seqs) {
sub_H <- subset(pf, V1 == i )
sub_D <- subset(pd, V1 == i )

# print(all.equal(sub_H,sub_D))
df <- resdf %>% filter(full.1 == i)
x <- unique(sort(df$taxa))

y <- gsub(x=i, pattern = "_.*", replacement="")
p <- ggplot() +
  geom_line(data = sub_H, aes(x=V2, y=log10(V3)), color = "#44CCD0")+
  geom_line(data = sub_D, aes(x=V2, y=log10(V3)), color="#FA918B")+
  scale_x_continuous(name="Base Position") +
  scale_y_continuous(name="log10 coverage")+
  ggtitle(x)+
  geom_vline(xintercept = sub_H$newA1, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newA2, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newB1, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newB2, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newC1, color = "black", size=.2)+
  geom_vline(xintercept = sub_H$newC2, color = "black", size=.2)+
  theme_minimal()
pdf(paste0(y, "_HI_cov",".pdf"))
print(p)
dev.off()
}

system("ls | grep _HI_cov | while read line; do ~/.iterm2/imgcat $line;done")
ureABC_HI_name <- names(which(table(sig_df$seqs) >= 3))
```