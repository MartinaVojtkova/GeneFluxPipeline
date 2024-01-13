library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Three arguments required: 1) Path to presence_absence_matrix.tsv, 2) path to the output directory, 3) slice size for gene selction", call. = FALSE)
}

PA_matrix_file <- args[1]
output_path <- args[2]
slice_size <- as.numeric(args[3])



####### FUNCTIONS #########

#calculate median of the samples where gene is present
cluster_median <- function(x) {
  gene_present <- x[x != 0]  
  if (length(gene_present) == 0) {
    return(0)  
  } else {
    return(median(gene_present))  
  }
}


##### CODE ######


#read presence_absence_matrix
presence_absence_matrix <- read.table(PA_matrix_file, 
                                      header = TRUE, 
                                      sep = "\t", 
                                      row.names = 1)
#convert to matrix 
presence_absence_matrix <- as.matrix(presence_absence_matrix)

gc()

#filter out genes present in less than 10 genomes 
presence_absence_matrix <- presence_absence_matrix[rowSums(presence_absence_matrix != 0) >= 10, ]


#Calculate median occurence
gene_median_occurrence <- data.frame(
  Gene = rownames(presence_absence_matrix),
  MedianOccurrence = apply(presence_absence_matrix, 
                           1, 
                           FUN = cluster_median)
)

#Select gene with median == 1 (present in each genome once)
genes_median_1 <- gene_median_occurrence %>% 
  filter(MedianOccurrence == 1.00) %>%
  select(Gene)

#Take every Xth gene
max_rows <- nrow(genes_median_1)
gene_list <- genes_median_1 %>%
  slice(seq(1, min(max_rows, 
                   n()), 
            by = slice_size)) 

#Create gene list 
gene_list %>%
  select(Gene) %>%
  write_tsv(file = file.path(output_path, "gene_list.txt"), 
            col_names = FALSE)