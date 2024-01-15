library(tidyverse)
library(vegan)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Three arguments required: 1) dRep sim matrix, 2) ANI matrix, 3) output file path", call. = FALSE)
}

sim_matrix_file <- args[1]
ani_file <- args[2]
output_file <- args[3]

gene <- str_extract(ani_file, "(?<=fastani_results\\/)[^\\/]+")


load_distance_matrix <- function(file_name) {
  lines <- readLines(file_name)
  len <- as.numeric(lines[1])+1

  if (len < 10) {
    return(NULL)  # Return NULL if len is less than 10
  }

  matrix_data <- read.table(text = lines[-1], fill = TRUE, col.names = 1:len, row.names = 1)
  seqNames <- make.unique(rownames(matrix_data))
  rownames(matrix_data) <- seqNames
  colnames(matrix_data) <- seqNames

  distance_matrix <- as.matrix(matrix_data)
  distance_matrix[upper.tri(distance_matrix)] <- t(distance_matrix)[upper.tri(distance_matrix)]
  diag(distance_matrix) <- 0

  return(distance_matrix)
}



# Load matrices
sim_genome_matrix <- read.table(sim_matrix_file)
sim_genome_matrix <- as.matrix(sim_genome_matrix)
ani_matrix <- load_distance_matrix(ani_file)

regex_pattern <- paste0("results/fastani_input/sequences/", gene, "/(.*)\\.fasta")

# Now use this pattern in gsub
colnames(ani_matrix) <- gsub(regex_pattern, "\\1", colnames(ani_matrix))
rownames(ani_matrix) <- gsub(regex_pattern, "\\1", rownames(ani_matrix))

ani_matrix <- ani_matrix/100

# Check if kma_dist_matrix is NULL
if (is.null(ani_matrix)) {
  result_df <- data.frame(gene_name = gene, statistic = NA, signif = NA)
  write.table(result_df, output_file, row.names = FALSE, col.names = FALSE) 
  quit(save = "no", status = 0) 
} else {
  # Reduce Genome matrix
  gene_presence <- rownames(ani_matrix)
  reduced_genome <- sim_genome_matrix[gene_presence, gene_presence]

  mantel_results <- mantel(ani_matrix,reduced_genome, method="pearson", permutations=999, strata = NULL,
                na.rm = TRUE, parallel = getOption("mc.cores"))

  result_df <- data.frame(
    gene_name = gene,
    statistic = mantel_results$statistic,
    signif =  mantel_results$signif
  )

  write.table(result_df, output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
}