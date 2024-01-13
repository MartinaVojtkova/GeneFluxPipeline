library(tidyverse)
library(vegan)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("4 arguments required: 1) dRep results, 2) input list metadata, 3) distance matrix output, 4) similiarity matrix output", call. = FALSE)
}

dRep_file <- args[1]
input_file <- args[2]
dist_output <- args[3]
sim_output <- args[4]

samples <- read.table(input_file)
samples$V2 <- sapply(strsplit(samples$V2, "/"), function(x) tail(x, 1))
dRep_results <- read_csv(dRep_file)

#match sample names to input 
result_df <- dRep_results %>%
  mutate(
    genome1 = samples$V1[match(genome1, samples$V2)],
    genome2 = samples$V1[match(genome2, samples$V2)]
  )

#make pairwise matrices 
dist_matrix <- acast(result_df, genome1 ~ genome2, value.var = "dist", fill = 0)
similarity_matrix <- acast(result_df, genome1 ~ genome2, value.var = "similarity", fill = 0)

#write files
write.table(dist_matrix, file = dist_output, sep = "\t", quote = FALSE, col.names = NA)
write.table(similarity_matrix, file = sim_output, sep = "\t", quote = FALSE, col.names = NA)