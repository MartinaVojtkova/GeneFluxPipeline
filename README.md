# GeneFluxPipeline

## Introduction 
This pipeline was designed to analyse the gene synteny and gene flux of core and accessory genes in bacteria. The input genomes are automatically annotated. Next, pan-genome analysis is performed, which serves as a base for gene selection. Selected sample of core and accessory genes is then subject to gene synteny analysis. The k-mer hamming distance, flank identity index  and average nucleotide identity matrices are calculated. Whole genome distances/similarities are determined. Distance differences between whole-genomes, target genes and their flanking environments are then analyzed using Procrustes Analysis. 

![Pipeline (2)](https://github.com/MartinaVojtkova/GeneFluxPipeline/assets/101507399/f6bffe38-e6a9-4289-ad55-659f794bc7f4)

## Input
A tab-separated file that links the sample names to their coresponding genome data file paths needs to be provided as input. Template of the file can be found in _Input/input_template.tsv_.
The pipeline's parameters can be changed in the config.yaml file. 

### Pipeline Parameters (config.yaml)
- **input_list:** List of the sample names and file paths. 
- **flank_length:** The length of the upstream and downstream flanking sequence of the target genes. (Default is 10 000)
- **drep_pc:** dRep primary clustering threshold for MASH (Default is 0.99) 
- **drep_sc:** dRep secondary clustering threshold (Default is 0.95) 
- **genus:**  Genus for bacterial genome Annotation
- **slice_size:** Portion of all genes to analyse (3 = one third of genes)
- **k-mer_size:** K-mer lenght for KMA distance calculation

## How to Run 
This is a snakemake pipeline. To execute Snakemake and Miniconda have to be installed. The default number of cores is 40. This will be downscaled if less is available. 
Run the following command: 
```
snakemake --use-conda --cores {CPUs}
```
## Output 
The pipeline provides genome annotation files from Prokka and the output of pan-genome analysis by Roary (excluding the temporary files)
The main output of the pipeline is gene_metrics.tsv file. 

**gene_metrics.tsv** contains: 
- Gene: Gene cluster name as provided by Roary 
- Kmer_distance_median: Per-gene median value of the flanking region k-mer distance matrix. (K-mer Hamming distance used as measure)
- ANI_median: Per-gene median value of the flanking region average nucleotide identity matrix
- Flank_index_median: Per-gene median value of the flanking region flank identity index matrix. (Measure of flank annotation similiarity)

## This Pipeline Uses 
