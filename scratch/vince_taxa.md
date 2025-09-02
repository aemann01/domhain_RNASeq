# 1. Get files from hillary
```sh
scp suzanne@hillary.clemson.edu:/home/allie/long_oral_microbiome/01-read_processing/rep_set.fa ~/rna_dohmain/scratch
scp suzanne@hillary.clemson.edu:/home/allie/long_oral_microbiome/01-read_processing/sequence_table.merged.txt ~/rna_dohmain/scratch
scp suzanne@hillary.clemson.edu:/home/allie/long_oral_microbiome/map_domhain_long.txt ~/rna_dohmain/scratch
```
# 2. Reassign taxonomy using homd database
```sh
#taxonomic assignment
kraken2 --db ../rpoc/database/kraken_homd \
	--threads 6 \
	--use-names \
	--output rep_set.kraken.out rep_set.fa \
	--unclassified-out rep_set.unclassified.out --confidence 0.01
# 	Loading database information... done.
# 38252 sequences (18.29 Mbp) processed in 0.406s (5654.0 Kseq/m, 2703.66 Mbp/m).
#   35684 sequences classified (93.29%)
#   2568 sequences unclassified (6.71%)
```
# 3. Make taxonomy file
```sh
awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > taxids
sed "s/\t//g" ../rpoc/database/rankedlineage.dmp > rankedlineage.dmp
sort -t "|" -k 1b,1 rankedlineage.dmp > rankedlineage_sorted
cat rankedlineage_sorted | sed 's/|\{2,\}/|/g' > rankedlineage_clean
# add unclassified to taxonomy file
sed -i '1 i\0|unclassified|' rankedlineage_clean
sed 's/|/\t/' rankedlineage_clean | sed 's/ /_/g' >rankedlineage_clean2
python3 ../rpoc/nt/lineages.py #output is lineage
awk -F"\t" '{print $2}' lineage | awk -F\| '{s=$NF;for(i=NF-1;i>=1;i--)s=s FS $i;print s}' | sed 's/^|//' | sed 's/ /_/g' | sed 's/|/;/g' > taxonomy
# merge assemebleies ids and taxonomy
awk '{print $2}' rep_set.kraken.out > asvids
paste asvids taxonomy > taxonomy.txt
#filter out what was only assigned at phylum level
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq ./rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python ~/bin/fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt
```
# 4. Add taxonomy onto Sequence table
transpose
```py
import pandas as pd
# Read the file into a DataFrame
df = pd.read_csv("sequence_table.merged.txt", sep='\t', index_col=0)  # Assuming SampleID is the first column
transposed_df = df.T
transposed_df.to_csv("trans.sequence_table.txt", sep='\t', index=True)
```
```sh
parallel -a <(awk '{print $1}' taxonomy_bac.txt) -j 7 -k "grep -wm 1 '{}' trans.sequence_table.txt " > sub.sequence_table.txt
paste -d "\t" <(awk '{print $8}' taxonomy_bac.txt) sub.sequence_table.txt > temp
cat <(head -n 1 trans.sequence_table.txt | sed 's/DM/species\tASV\tDM/') temp > temp2
mv temp2 sub.sequence_table.txt
sed 's/species\t//' sub.sequence_table.txt > test
```
transpose again
```py
import pandas as pd
# Read the file into a DataFrame
df = pd.read_csv("test", sep='\t', index_col=0)  # Assuming SampleID is the first column
df.index.name = "species"  # Remove the name of the index
transposed_df = df.T
transposed_df.to_csv("taxa.sequence_table.txt", sep='\t', index=True, header=True)
```
Get health cateogory
```sh
sed -i 's/Lachnospiraceae/species\tLachnospiraceae/' taxa.sequence_table.txt
grep DM00 taxa.sequence_table.txt > samples.sequence_table.txt
parallel -a <(awk '{print $1}' samples.sequence_table.txt) -j 7 -k "grep -wm 1 '{}' map_domhain_long.txt" | awk '{print $18}' > tooth_health
paste -d "\t" tooth_health samples.sequence_table.txt > tooth.sequence_table.txt
head -n 2 taxa.sequence_table.txt > species_asvs
sed -i 's/species/health\tspecies/' species_asvs
sed -i 's/ASV/health\tASV/' species_asvs
cat species_asvs tooth.sequence_table.txt > total.sequence_table.txt
```
transpose again
```py
import pandas as pd
# Read the file into a DataFrame
df = pd.read_csv("total.sequence_table.txt", sep='\t', index_col=0)  # Assuming SampleID is the first column
transposed_df = df.T
transposed_df.to_csv("taxa.sequence_table.txt", sep='\t', index=True, header=True)
```
```sh
sed -i 's/health/health\thealth/' taxa.sequence_table.txt
```
# 4. Make one collpased in R
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
setwd("~/rna_dohmain/scratch")
seqtab <- t(read.table("./sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("./taxonomy_bac.txt", header=F, row.names=1, sep="\t")
# tree <- read.tree("RAxML_bestTree.ref.tre")
# tree.root <- midpoint.root(tree)
map <- read.table("./map_domhain_long.txt", sep="\t", header=T, row.names=1)
notinmeta <- setdiff(colnames(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), colnames(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw

ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
ps.dat
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
long_df <- reshape2::dcast(data, Sample ~ V8, value.var = "Abundance", fun.aggregate = sum)
long_df$Sample <- as.character(long_df$Sample)
map$Sample <- row.names(map)
final_df <- merge(long_df, map[, c("Sample", "tooth_health")], by="Sample", all.x=TRUE)
last_col <- names(final_df)[ncol(final_df)]
# Reorder the dataframe with the last column as the 2nd column
final_df <- final_df[, c("Sample", last_col, setdiff(names(final_df), c("Sample", last_col)))]

widedf <- t(final_df)
write.table(widedf, "transposed_data.txt", sep = "\t", quote = FALSE, row.names = TRUE)
```
```sh
sed -i '1d' transposed_data.txt
mv transposed_data.txt species.sequence_table.txt
```
# Make taxonomy for US on hillary
```sh
cd /home/allie/ads_plaque/05-rpoC_processing/ ~
cp sequence_table.merged.txt ~/rna_dohmain/scratch
cp taxonomy_bac.txt ~/rna_dohmain/scratch
cp map.txt ~/rna_dohmain/scratch
```
## 1.  Add taxonomy onto Sequence table
transpose
```py
import pandas as pd
# Read the file into a DataFrame
df = pd.read_csv("sequence_table.merged.txt", sep='\t', index_col=0)  # Assuming SampleID is the first column
transposed_df = df.T
transposed_df.to_csv("trans.sequence_table.txt", sep='\t', index=True)
```
```sh
parallel -a <(awk '{print $1}' taxonomy_bac.txt) -j 50 -k "grep -wm 1 '{}' trans.sequence_table.txt " > sub.sequence_table.txt
paste -d "\t" <(awk '{print $10}' taxonomy_bac.txt) sub.sequence_table.txt > temp
cat <(head -n 1 trans.sequence_table.txt | sed 's/DM/species\tASV\tDM/') temp > temp2
mv temp2 sub.sequence_table.txt
sed 's/species\t//' sub.sequence_table.txt > test
```
transpose again
```py
import pandas as pd
# Read the file into a DataFrame
df = pd.read_csv("test", sep='\t', index_col=0)  # Assuming SampleID is the first column
df.index.name = "species"  # Remove the name of the index
transposed_df = df.T
transposed_df.to_csv("taxa.sequence_table.txt", sep='\t', index=True, header=True)
```
Get health cateogory
```sh
sed -i 's/Prevotella_pallens/species\tPrevotella_pallens/' taxa.sequence_table.txt
grep UF taxa.sequence_table.txt > samples.sequence_table.txt
awk '{print $1}' samples.sequence_table.txt | sed 's/.*P/P/' | sed 's/R//' > tooth_health
paste -d "\t" tooth_health samples.sequence_table.txt > tooth.sequence_table.txt
head -n 2 taxa.sequence_table.txt > species_asvs
sed -i 's/species/health\tspecies/' species_asvs
sed -i 's/ASV/health\tASV/' species_asvs
cat species_asvs tooth.sequence_table.txt > total.sequence_table.txt
```
transpose again
```py
import pandas as pd
# Read the file into a DataFrame
df = pd.read_csv("total.sequence_table.txt", sep='\t', index_col=0)  # Assuming SampleID is the first column
transposed_df = df.T
transposed_df.to_csv("taxa.sequence_table.txt", sep='\t', index=True, header=True)
```
```sh
sed -i 's/health/health\thealth/' taxa.sequence_table.txt
```
## 2. Get sequence table by species
Make taxonomy file with all samples
```sh
paste -d "\t" <(awk '{print $1}' samples.sequence_table.txt) tooth_health > map.txt
sed -i '1 i\sample\ttooth_health' map.txt
```
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
setwd("~/rna_dohmain/scratch")
seqtab <- t(read.table("./sequence_table.merged.txt", header=T, row.names=1))
tax <- read.table("./taxonomy_bac.txt", header=F, row.names=1, sep="\t")
# tree <- read.tree("RAxML_bestTree.ref.tre")
# tree.root <- midpoint.root(tree)
map <- read.table("./map.txt", sep="\t", header=T, row.names=1)
notinmeta <- setdiff(colnames(seqtab), row.names(map))
notinraw <- setdiff(row.names(map), colnames(seqtab))
print("Samples found in ASV table but not in metadata:")
notinmeta
print("Samples found in metadata but not in sequencing table:")
notinraw

ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=T), sample_data(map), tax_table(as.matrix(tax)))
ps.dat
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[8])
data <- psmelt(glom) # create dataframe from phyloseq object
long_df <- reshape2::dcast(data, Sample ~ V9, value.var = "Abundance", fun.aggregate = sum)
long_df$Sample <- as.character(long_df$Sample)
map$Sample <- row.names(map)
final_df <- merge(long_df, map[, c("Sample", "tooth_health")], by="Sample", all.x=TRUE)
last_col <- names(final_df)[ncol(final_df)]
# Reorder the dataframe with the last column as the 2nd column
final_df <- final_df[, c("Sample", last_col, setdiff(names(final_df), c("Sample", last_col)))]

widedf <- t(final_df)
write.table(widedf, "transposed_data.txt", sep = "\t", quote = FALSE, row.names = TRUE)
```
```sh
sed -i '1d' transposed_data.txt
mv transposed_data.txt species.sequence_table.txt