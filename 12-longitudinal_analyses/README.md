

```bash
cd ~/domhain_RNAseq/09-longitudinal_analyses
conda activate 2024-HIV_RNASeq



```R
metadata <- read.table("/home/allie/domhain_RNAseq/map.txt", header=T, sep="\t")
# remove dashes from health categories or it will mess up downstream processing
metadata$aliquot_type <- sub("-", "", metadata$aliquot_type)
row.names(metadata) <- metadata$sample_id
# read in gene counts file
genecounts <- read.table("/home/allie/domhain_RNAseq/03-star_map/02-HOMD_map/featurecounts/read_counts.txt", header=T, sep="\t", row.names=1)
# get rid of weird empty column in genecounts
genecounts <- genecounts[1:(length(genecounts)-1)]
# fix sample names in gene counts so they match the metadata
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "\\.", replacement = "-")  
# reorder columns by metadata 
metadata <- metadata[order(colnames(genecounts)),]
# colnames(genecounts)
# rownames(metadata)
# check to make sure that sample ids match between gene counts and metadata
table(colnames(genecounts)==metadata$sample_id) # should return all true