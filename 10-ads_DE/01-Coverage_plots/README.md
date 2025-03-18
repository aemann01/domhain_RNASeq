# Coverage plots of highly expressed species

I want to see if these are similar to what I was seeing with the UF data (i.e., steep coverage dropoff at arcC)

First pull the operons that I'm most interested in

```bash
# from all H v D deseq data, pull genes that had a fold change of at least 2 and an adjusted p value of 0.05 or less. For each genome at least three genes must be upregulated
cd /home/allie/domhain_RNAseq/07-ads_expression
mkdir coverage_plots && cd coverage_plots
awk -F"\t" '$2 >= 2 || $2 <= -2' ../deseq_results_ADS-HvD.txt | awk -F"\t" '$5 <= 0.05' | awk -F"\t" '{print $1}' | sed 's/_.*//' | sort | uniq -c | sed 's/^[ \t]*//g' | awk -F" " '$1 >= 3' | awk '{print $2}' > genomes_of_interest.ids
# filter out these operons from the sliced operon file
cat genomes_of_interest.ids| while read line; do grep $line -m1 ../ads_operons.fna -A 1; done > genomes_of_interest.fa

# First need to build a database from the HOMD operons
conda install bioconda::bowtie2
bowtie2-build genomes_of_interest.fa genomes_of_interest.db
# move to folder containing filtered data
cd ~/domhain_RNAseq/01-processing/filtered
# map to each genome
ls *.1.fastq.gz | sed 's/.1.fastq.gz//' | parallel -j 50 --gnu 'bowtie2 -x ~/domhain_RNAseq/07-ads_expression/coverage_plots/genomes_of_interest.db -1 {}.1.fastq.gz -2 {}.2.fastq.gz --end-to-end  --qc-filter --no-unal -t -S ~/domhain_RNAseq/07-ads_expression/coverage_plots/{}.sam 2>~/domhain_RNAseq/07-ads_expression/coverage_plots/{}.out 1>~/domhain_RNAseq/07-ads_expression/coverage_plots/{}.err'

cd ~/domhain_RNAseq/07-ads_expression/coverage_plots
# convert to bam and sort by coordinate
ls *sam | sed 's/.sam//' | while read line; do samtools sort $line.sam -o $line.sort.bam
; done
# split bam files by reference
ls *sort.bam | sed 's/.sort.bam//' | while read line; do bamtools split -in $line.sort.bam -reference; done
# clean up file names
for file in *sort*.bam; 
	do dest="${file//[[:space:]]/.}" && mv -i "$file" "${dest//[^[:alnum:]._-]/}"; 
done
# merge based on reference
for ref in $(ls *SEQF*bam | rev | cut -f1 -d_ | rev | sort -u)
	do cat *$ref > $ref
done
# get depth by position for each reference
ls SEQ* | sed 's/.bam//' | while read line; do samtools depth $line.bam > $line.depth; d
one		
```

Plot in R

```R
# test with sanguinis first
seqf1070 <- read.table("SEQF1070.1CP000387.1718765-722174.depth", header=F, sep="\t")
pdf("seqf1070_coverage.pdf")
plot(seqf1070$V2, log10(seqf1070$V3), type="l", xlab="base position", ylab="log10 coverage", lwd=3)
dev.off()
# this is really wobbly, try with other species to see if we have any with good coverage across the operon
filelist <- list.files(pattern="*depth")
for(i in 1:length(filelist)){
	dat <- read.table(filelist[i], header=F, sep="\t")
	pdf(paste(filelist[i], "_cov.plot.pdf", sep=""))
	plot(dat$V2, log10(dat$V3), type="l", xlab="base position", ylab="log10 coverage",lwd=3, col="#FA918B")
dev.off()
}
# none of these have very good coverage - I wonder if we just don't have the sequencing depth.
```