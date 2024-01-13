library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Three arguments required: 1) Path to the input list and metadata, 2) path to the Roary clustered_proteins file, 3) path to the output directory", call. = FALSE)
}

input_list_metadata <- args[1]
cluster_file <- args[2]
output_path <- args[3]
output_file <- paste0(output_path,"presence_absence_matrix.tsv")

###### CODE #####

#Load list of assemblies
assembly_names <- read.table(input_list_metadata, 
                             header = FALSE, 
                             sep = "\t")
assembly_names <- assembly_names$V1


#open files for r/w
infile <- file(cluster_file, 
               open = "r")
outfile <- file(output_file, 
                open = "w")


#write the assembly names as header 
writeLines(paste(c("Gene", 
                   assembly_names), 
                 collapse = "\t"), 
           outfile)

#set up a counter 
assembly_counter <- setNames(rep(0, length(assembly_names)), 
                          assembly_names)


while(length(line <- readLines(infile, n = 1, warn = FALSE)) > 0) {
  
  parts <- strsplit(line, ": ", fixed = TRUE)[[1]]    #split gene_name and genomes 
  gene <- parts[1]                                    #save gene_name
  genomes <- strsplit(parts[2], "\t")[[1]]            #save genomes
  genomes <- gsub("_\\d+$", "", genomes)              #remove trailing indexes
  
  #Count occurrences
  counts <- table(factor(genomes, 
                         levels = names(assembly_counter)))
  
  #Update the dictionary with counts for present assemblies
  assembly_counter[names(counts)] <- counts
  
  #writte counts to outfile
  output_line <- paste(gene, 
                       paste(assembly_counter, 
                                   collapse = "\t"), 
                       sep = "\t")
  
  writeLines(output_line, outfile)
  
  #reset counter to zero
  assembly_counter <- setNames(rep(0, length(assembly_names)), 
                               assembly_names)
}


close(infile)
close(outfile)
