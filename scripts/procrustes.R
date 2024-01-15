library(tidyverse)
library(vegan)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Three arguments required: 1) dRep distance matrix, 2) k-mer distance matrix, 3) output file path", call. = FALSE)
}

dist_matrix_file <- args[1]
kma_dist_file <- args[2]
output_file <- args[3]

gene <- str_extract(kma_dist_file, "(?<=kma_dist\\/)[^\\/]+")


load_distance_matrix <- function(file_name) {
    lines <- readLines(file_name)
    len <- as.numeric(lines[1]) +1
    matrix_data <- read.table(text = lines[-1], fill = TRUE, col.names = 1:len, row.names = 1)

    seqNames <- make.unique(rownames(matrix_data))
    rownames(matrix_data) <- seqNames
    colnames(matrix_data) <- seqNames

    distance_matrix <- as.matrix(matrix_data)

    distance_matrix[upper.tri(distance_matrix)] <- t(distance_matrix)[upper.tri(distance_matrix)]
    diag(distance_matrix) <- 0

    return(distance_matrix)
}

normalize_matrix <- function(matrix) {
  min_val <- min(matrix)
  max_val <- max(matrix)
  (matrix - min_val) / (max_val - min_val)
}


dist_genome_matrix <- read.table(dist_matrix_file)
dist_genome_matrix <- as.matrix(dist_genome_matrix)
kma_dist_matrix <- load_distance_matrix(kma_dist_file)


gene_presence <- rownames(kma_dist_matrix)
reduced_genome <- dist_genome_matrix[gene_presence, gene_presence]

#ordination
kma_norm <- normalize_matrix(kma_dist_matrix)
ord_genome <- cmdscale(reduced_genome)
ord_kma <- cmdscale(kma_norm)

 procur_result <- protest(ord_genome, ord_kma)

result_df <- data.frame(
  gene_name = gene,
  cor = procur_result$t0,
  signif =  procur_result$signif
  )

write.table(result_df, output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
