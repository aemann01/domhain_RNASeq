# 1. Redo feeaturecounts for HOMD to get the read IDs that aligned
```sh
cd ~/rna_dohmain/03-star_map/02-HOMD_map
ls *Aligned.out.bam | sed 's/Aligned.out.bam//' | while read line; do featureCounts -f -p -C -B -M -a ../homd_db/ALL_genomes.gtf -o featurecounts/$line.out -R BAM -T 60 $line\Aligned.out.bam -t CDS -g transcript_id; done
```
# 2. Make list of reads that aligned to SEQF6415.1
```sh
cd featurecounts
ls *.homdAligned.out.bam.featureCounts.bam | sed 's/.homdAligned.out.bam.featureCounts.bam//' | parallel -j 23 'samtools sort -@ 8 -o {1}.sorted.bam {1}.homdAligned.out.bam.featureCounts.bam' # sort the bam files
ls *.sorted.bam | while read line; do samtools index -@ 150 $line; done # index the files
# get arcA (cords find using ALL_genomes.gtf)
parallel -j 93 'samtools view -h {} "SEQF6415.1|CP046524.1:1855788-1857017" > {.}_arcA.sam' ::: *.sorted.bam
# get arcB
parallel -j 93 'samtools view -h {} "SEQF6415.1|CP046524.1:1857074-1858090" > {.}_arcB.sam' ::: *.sorted.bam
# get arcC
parallel -j 93 'samtools view -h {} "SEQF6415.1|CP046524.1:1858384-1859331" > {.}_arcC.sam' ::: *.sorted.bam

# get the read IDs for arcABC
parallel -j 93 'base={.};grep SEQF6415 {} | grep XT | awk "{print \$1}" | sort | uniq > ${base}.ids' ::: *.sorted_arcA.sam
parallel -j 93 'base={.};grep SEQF6415 {} | grep XT | awk "{print \$1}" | sort | uniq > ${base}.ids' ::: *.sorted_arcB.sam
parallel -j 93 'base={.};grep SEQF6415 {} | grep XT | awk "{print \$1}" | sort | uniq > ${base}.ids' ::: *.sorted_arcC.sam
```
# 3. Redo feature counts for denovo to get the read IDs that aligned
```sh
cd ~/rna_dohmain/12-denovo_analyses
ls *Aligned.out.bam | sed 's/Aligned.out.bam//' | parallel -j 10 '
  sample={};
  featureCounts -f -p -C -B -M \
    -a ./arcABC_assemblies.gtf \
    -o featurecounts/${sample}.out \
    -R BAM -T 19 \
    ${sample}Aligned.out.bam \
    -t CDS -g transcript_id
'
```
# 4. Sort and get the IDS for the denovo
```sh
cd featurecounts
ls *.denovoAligned.out.bam.featureCounts.bam | sed 's/.denovoAligned.out.bam.featureCounts.bam//' | parallel -j 23 'samtools sort -@ 8 -o {1}.sorted.bam {1}.denovoAligned.out.bam.featureCounts.bam' # sort the bam 
ls *.sorted.bam | while read line; do samtools index -@ 150 $line; done # index the files

# convert bam to sam
# get arcA
grep arcA ../arcABC_assemblies.gtf | awk '$3=="transcript" {
  chrom=$1; start=$4-1; end=$5;
  name="arcA_" NR;
  match($0, /gene_name "([^"]+)"/, arr);
  if (arr[1]) name=arr[1];
  print chrom"\t"start"\t"end"\t"name
}' > arcA_regions.bed
ls *.sorted.bam | parallel -j 93 '
  bedtools intersect -abam {} -b arcA_regions.bed > {.}_arcA.bam && \
  samtools index {.}_arcA.bam && \
  samtools view -h {.}_arcA.bam > {.}_arcA.sam
'
# get arcB
grep arcB ../arcABC_assemblies.gtf | awk '$3=="transcript" {
  chrom=$1; start=$4-1; end=$5;
  name="arcB_" NR;
  match($0, /gene_name "([^"]+)"/, arr);
  if (arr[1]) name=arr[1];
  print chrom"\t"start"\t"end"\t"name
}' > arcB_regions.bed
ls *.sorted.bam | parallel -j 93 '
  bedtools intersect -abam {} -b arcB_regions.bed > {.}_arcB.bam && \
  samtools index {.}_arcB.bam && \
  samtools view -h {.}_arcB.bam > {.}_arcB.sam
'
# get arcC
grep arcC ../arcABC_assemblies.gtf | awk '$3=="transcript" {
  chrom=$1; start=$4-1; end=$5;
  name="arcC_" NR;
  match($0, /gene_name "([^"]+)"/, arr);
  if (arr[1]) name=arr[1];
  print chrom"\t"start"\t"end"\t"name
}' > arcC_regions.bed
ls *.sorted.bam | parallel -j 93 '
  bedtools intersect -abam {} -b arcC_regions.bed > {.}_arcC.bam && \
  samtools index {.}_arcC.bam && \
  samtools view -h {.}_arcC.bam > {.}_arcC.sam
'
# now see which read IDs from that sample map to a denovo sequence for arcABC
#arcA
ls *.sorted_arcA.sam | parallel -j 93 '
  sample={.}
  grep -Fwf ../../03-star_map/02-HOMD_map/featurecounts/${sample}.ids {} | grep XT | awk "{print \$3}" | sort | uniq -c > ${sample}.arcA.count
  grep -Fwf ../../03-star_map/02-HOMD_map/featurecounts/${sample}.ids {} | grep XT | awk "{print \$1}" | sort | uniq -c | awk "{print \$2}" > ${sample}.ids.found
'
#arcB
ls *.sorted_arcB.sam | parallel -j 93 '
  sample={.}
  grep -Fwf ../../03-star_map/02-HOMD_map/featurecounts/${sample}.ids {} | grep XT | awk "{print \$3}" | sort | uniq -c > ${sample}.arcB.count
  grep -Fwf ../../03-star_map/02-HOMD_map/featurecounts/${sample}.ids {} | grep XT | awk "{print \$1}" | sort | uniq -c | awk "{print \$2}" > ${sample}.ids.found
'
#arcC
ls *.sorted_arcC.sam | parallel -j 93 '
  sample={.}
  grep -Fwf ../../03-star_map/02-HOMD_map/featurecounts/${sample}.ids {} | grep XT | awk "{print \$3}" | sort | uniq -c > ${sample}.arcC.count
  grep -Fwf ../../03-star_map/02-HOMD_map/featurecounts/${sample}.ids {} | grep XT | awk "{print \$1}" | sort | uniq -c | awk "{print \$2}" > ${sample}.ids.found
'

# find the missing ids 
#arcA
ls *sorted_arcA.ids.found | sed 's/.sorted_arcA.ids.found//' | while read line; do cat $line.sorted_arcA.ids.found ../../03-star_map/02-HOMD_map/featurecounts/$line.sorted_arcA.ids | sort | uniq -u > $line.sorted_arcA.ids.not_found; done
#arcB
ls *sorted_arcB.ids.found | sed 's/.sorted_arcB.ids.found//' | while read line; do cat $line.sorted_arcB.ids.found ../../03-star_map/02-HOMD_map/featurecounts/$line.sorted_arcB.ids | sort | uniq -u > $line.sorted_arcB.ids.not_found; done
#arcC
ls *sorted_arcC.ids.found | sed 's/.sorted_arcC.ids.found//' | while read line; do cat $line.sorted_arcC.ids.found ../../03-star_map/02-HOMD_map/featurecounts/$line.sorted_arcC.ids | sort | uniq -u > $line.sorted_arcC.ids.not_found; done

# find the amount misisng
ls *.arcA.count | parallel -j 93 "
  sample={= s/\\.arcA\\.count\$// =}
  total_ids=\$(wc -l < ../../03-star_map/02-HOMD_map/featurecounts/\${sample}.ids)
  found_ids=\$(wc -l < \${sample}.ids.found)
  percent_found=\$(awk -v f=\"\$found_ids\" -v t=\"\$total_ids\" 'BEGIN { if(t>0) printf \"%.2f\", (f/t)*100; else print \"NA\" }')
  echo -e \"\${sample}\t\${found_ids}\t\${total_ids}\t\${percent_found}%\"
" > reads_found_arcA.tsv

ls *.arcB.count | parallel -j 93 "
  sample={= s/\\.arcB\\.count\$// =}
  total_ids=\$(wc -l < ../../03-star_map/02-HOMD_map/featurecounts/\${sample}.ids)
  found_ids=\$(wc -l < \${sample}.ids.found)
  percent_found=\$(awk -v f=\"\$found_ids\" -v t=\"\$total_ids\" 'BEGIN { if(t>0) printf \"%.2f\", (f/t)*100; else print \"NA\" }')
  echo -e \"\${sample}\t\${found_ids}\t\${total_ids}\t\${percent_found}%\"
" > reads_found_arcB.tsv

ls *.arcC.count | parallel -j 93 "
  sample={= s/\\.arcC\\.count\$// =}
  total_ids=\$(wc -l < ../../03-star_map/02-HOMD_map/featurecounts/\${sample}.ids)
  found_ids=\$(wc -l < \${sample}.ids.found)
  percent_found=\$(awk -v f=\"\$found_ids\" -v t=\"\$total_ids\" 'BEGIN { if(t>0) printf \"%.2f\", (f/t)*100; else print \"NA\" }')
  echo -e \"\${sample}\t\${found_ids}\t\${total_ids}\t\${percent_found}%\"
" > reads_found_arcC.tsv
```
# 5. Combine
```sh
# arcA
ls *.sorted_arcA.arcA.count | parallel -j 93 'sample={.}; awk -v s="$sample" '\''{print $2, $1, "arcA", s}'\'' {}  | sed "s/.sorted_arcA.arcA//"> ${sample}.arcA.counts'
#arcB
ls *.sorted_arcB.arcB.count | parallel -j 93 'sample={.}; awk -v s="$sample" '\''{print $2, $1, "arcB", s}'\'' {}  | sed "s/.sorted_arcB.arcB//"> ${sample}.arcB.counts'
#arcC
ls *.sorted_arcC.arcC.count | parallel -j 93 'sample={.}; awk -v s="$sample" '\''{print $2, $1, "arcC", s}'\'' {}  | sed "s/.sorted_arcC.arcC//"> ${sample}.arcC.counts'
# combine all
cat *counts > arcABC.counts
```
# 6. Make heatmap in R
```R
library(DESeq2, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)
library(EnhancedVolcano, warn.conflicts = F, quietly = T)
library(dplyr, warn.conflicts = F, quietly = T)
library(viridis, warn.conflicts = F, quietly = T)
library(phyloseq, warn.conflicts = F, quietly = T)
library(phytools, warn.conflicts = F, quietly = T)
library(phyloseq.extended, warn.conflicts = F, quietly = T)
library(tidyr)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(patchwork)
library("cowplot")
require(phylobase)
library(phangorn)
library(stringr)
library(reshape2)
setwd("~/rna_dohmain/12-denovo_analyses")

# read in the file
metadata <- read.table("~/rna_dohmain/homd_map/map.txt", header=T, sep="\t")
read_counts <- read.table("~/rna_dohmain/12-denovo_analyses/featurecounts/arcABC.counts", header=F, sep=" ")
colnames(read_counts) <- c("locus_tag", "count", "gene", "sample")
genecounts <- read.csv("~/rna_dohmain/07-ads_expression/arcGene_read_counts.cleaned.txt", header=T, sep="\t", row.names=1)
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-") 
# find which locus tags were arcABC for superstar
homd <- read.csv("~/rna_dohmain/homd_map/annotations.arc.txt", header=T, sep="\t", quote="")
sub_homd <- dplyr::filter(homd, SEQ_ID == "SEQF6415.1")

# make heatmap
full_grid <- expand.grid(
  sample = metadata$sample_id,
  gene = c("arcA", "arcB", "arcC")
)
read_counts_complete <- full_grid %>%
  left_join(read_counts, by = c("sample", "gene")) %>%
  mutate(count = ifelse(is.na(count), 0, count))

samples <- unique(read_counts_complete$sample)
genes   <- unique(read_counts_complete$gene)
loci    <- unique(read_counts_complete$locus_tag)
full_grid <- expand_grid(locus_tag = loci, gene = genes, sample = samples)
read_counts_complete2 <- full_grid %>%
  left_join(read_counts_complete, by = c("locus_tag", "gene", "sample")) %>%
  mutate(count = ifelse(is.na(count), 0, count))

# add red astricks 
read_counts_complete2 <- read_counts_complete2 %>%
  group_by(gene, sample) %>%
  mutate(is_top = count == max(count, na.rm = TRUE) & count > 0) %>%
  ungroup()


pdf("denovo_mapped_SEQF6415.1.pdf", width =40, height =15)
ggplot(read_counts_complete2, aes(x = sample, y = locus_tag, fill = log10(count))) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  geom_point(data = subset(read_counts_complete2, is_top),
             aes(x = sample, y = locus_tag),
             color = "red", shape = 8, size = 2) +  # asterisk
  theme_minimal() +
  labs(title = "Top Locus per Gene-Sample (Red Asterisk)",
       x = "Sample", y = "Locus Tag", fill = "log10(Count + 1)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~ gene, scales = "free_x")
dev.off()
system("~/.iterm2/imgcat denovo_mapped_SEQF6415.1.pdf")



# compare to the original genome
summed_counts <- read_counts %>%
  group_by(gene, sample) %>%
  summarise(count = sum(count), .groups = 'drop') %>%
  complete(sample, gene = c("arcA", "arcB", "arcC"), fill = list(count = 0))
summed_counts <- summed_counts %>%
  mutate(sample_id = str_replace_all(sample, "\\.", "-"))
full_grid <- expand.grid(
  sample = metadata$sample_id,
  gene = c("arcA", "arcB", "arcC")
)
summed_counts_complete <- full_grid %>%
  left_join(summed_counts, by = c("sample", "gene")) %>%
  mutate(count = ifelse(is.na(count), 0, count))
summed_counts_complete$locus_tag <- "all_denovo"
summed_counts_complete$sample_id <- NULL
subcount <- genecounts[row.names(genecounts) %in% sub_homd$locus_tag, ]
subcount_long <- subcount %>%
  tibble::rownames_to_column(var = "locus_tag") %>%
  pivot_longer(
    cols = -locus_tag,
    names_to = "sample",
    values_to = "count"
  )
subcount_long2 <- left_join(as.data.frame(subcount_long), as.data.frame(sub_homd), by = c("locus_tag"))
subcount_long2 <- subcount_long2 %>%
  select(1:4) %>% 
  mutate(gene = ifelse(gene == "arcC1", "arcC", gene) )
subcount_long2$locus_tag <- "SEQF6415.1"
combined_counts <- rbind(summed_counts_complete, subcount_long2)
counts_all <- left_join(as.data.frame(combined_counts), metadata, by = c("sample" = "sample_id"))
counts_all$hiv_status <- factor(counts_all$hiv_status, levels=c("HI", "HEU", "HUU"))

# mark if superstar is higher
all_denovo_counts <- counts_all %>%
  filter(locus_tag == "all_denovo") %>%
  select(sample, gene, all_denovo_count = count)
counts_all2 <- counts_all %>%
  left_join(all_denovo_counts, by = c("sample", "gene")) %>%
  mutate(is_higher_than_denovo = (count > all_denovo_count) & (locus_tag != "all_denovo"))

pdf("SEQF6415.1vALL_denovo.pdf", width = 20)
ggplot(counts_all2, aes(x = sample, y = locus_tag, fill = log10(count))) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  geom_point(data = subset(counts_all2, is_higher_than_denovo),
             aes(x = sample, y = locus_tag),
             color = "red", size = 2, shape = 8) +  # star marker
  theme_minimal() +
  labs(title = "Read Count for SEQF6415.1 in Denovo ",
       x = "Sample",
       y = "Reference v Denovo",
       fill = "log10(Count)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~ hiv_status + gene, scales = "free_x")
dev.off()
system("~/.iterm2/imgcat SEQF6415.1vALL_denovo.pdf")

# see if overall per gene if super star is higher
summed_by_gene_locus <- counts_all2 %>%
  group_by(gene, locus_tag) %>%
  summarise(total_count = sum(count), .groups = "drop")

head(summed_by_gene_locus)
```
