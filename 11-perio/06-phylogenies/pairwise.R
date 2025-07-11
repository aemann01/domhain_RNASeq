# Load the seqinr library
library(seqinr, warn.conflicts = F, quietly = T)
library(matrixStats, warn.conflicts = F, quietly = T)

#get list of files
working_dir <- getwd()
file_list <- list.files(path = working_dir, 
                        pattern = "\\.align\\.fa$", 
                        full.names = TRUE)
# make empty list
distance_matrices <- list()

# loop through each alignment
for (file in file_list) {
  alignment <- read.alignment(file = file, format = "fasta")
  # print(paste("Processing file:", file))
  distance_matrix <- dist.alignment(alignment, matrix = "identity")
  distance_matrices[[basename(file)]] <- distance_matrix
}

# flatten list
df_distances <- as.data.frame(do.call(cbind, distance_matrices))
core_distances <- colMedians(as.matrix(df_distances))
# Read the alignment file
paste("Median pairewise distance of all the clusters:", median(core_distances))