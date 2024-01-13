# Gene Flux Analysis Pipeline

## Introduction 
This pipeline was designed to analyse the gene synteny and gene flux of core and accessory genes in bacteria. The input genomes are automatically annotated. Next, pan-genome analysis is performed, which serves as a base for gene selection. Selected sample of core and accessory genes is then subject to gene synteny analysis. The k-mer hamming distance, flank identity index  and average nucleotide identity matrices are calculated. Whole genome distances/similarities are determined. Distance differences between whole-genomes, target genes and their flanking environments are then analyzed using Procrustes Analysis. 

![Pipeline (2)](https://github.com/MartinaVojtkova/GeneFluxPipeline/assets/101507399/f6bffe38-e6a9-4289-ad55-659f794bc7f4)

## Input
A tab-separated file that links the sample names to their coresponding genome data file paths needs to be provided as input. Template of the file can be found [here](templates/input_template.tsv).
The pipeline's parameters can be changed in the [config.yaml](config.yaml) file. 

### Pipeline Parameters (config.yaml)
- **input_list:** List of the sample names and file paths. 
- **flank_length:** The length of the upstream and downstream flanking sequence of the target genes. (Default is 10 000)
- **drep_pc:**  `dRep compare` primary clustering threshold for MASH (Default is 0.99) 
- **drep_sc:** `dRep compare` secondary clustering threshold (Default is 0.95) 
- **genus:**  Genus for bacterial genome Annotation
- **slice_size:** Portion of all genes to analyse (3 = one third of genes)
- **k-mer_size:** K-mer length for KMA distance calculation

To choose the correct primary and secondary clustering thresholds, please consult the [dRep docummentation](https://drep.readthedocs.io/en/latest/overview.html).

## Run on Computerome
This is a snakemake pipeline that was designed and tested on Computerome(The Danish National Supercomputer for Life Sciences). 
To execute **Snakemake (v 6.1.9)** and **Miniconda (v 4.11.0)** have to be installed. The default number of cores is 40. This will be downscaled if less is available. 
Load the required packages: 
```
module purge
module load tools
module load miniconda3/4.11.0  
module load snakemake/6.9.1
```
Run the pipeline: 
```
snakemake --use-conda --cores 40 --max-status-checks-per-second 0.01
```
## Output 
The pipeline provides genome annotation files from Prokka and the output of pan-genome analysis by Roary (excluding the temporary files)
The main output of the pipeline is gene_metrics.tsv file. See example output [here](templates/gene_metrics_example.tsv). 

`gene_metrics.tsv` contains: 
- **Gene:** Gene cluster name as provided by Roary 
- **Kmer_distance_median:** Per-gene median value of the flanking region k-mer distance matrix. (K-mer Hamming distance used as measure)
- **ANI_median:** Per-gene median value of the flanking region average nucleotide identity matrix
- **Flank_index_median:** Per-gene median value of the flanking region flank identity index matrix. (Measure of flank annotation similiarity)
- **Proc_corelation:** The correlation between the flanking region k-mer distances and whole genome distances. (Symmetric Procrustes Analysis)
- **Proc_significance:** The significance value for Procrustes correlation. 

The complete distance/similiarity matrices are provided as supplementary output for each gene.  

## This Pipeline Uses 
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.09-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![dRep](https://img.shields.io/badge/dRep-≥3.4.5-brightgreen.svg?style=flat)](https://drep.readthedocs.io/en/latest/)
[![Prokka](https://img.shields.io/badge/Prokka-≥1.14.6-brightgreen.svg?style=flat)](https://github.com/tseemann/prokka)
[![Roary](https://img.shields.io/badge/Roary-≥3.13.0-brightgreen.svg?style=flat)](https://sanger-pathogens.github.io/Roary/)

[![Bedtools](https://img.shields.io/badge/Bedtools-≥2.27.1-brightgreen.svg?style=flat)](https://bedtools.readthedocs.io/en/latest/)
[![seqkit](https://img.shields.io/badge/seqkit-≥2.1.0-brightgreen.svg?style=flat)](https://bioinf.shenwei.me/seqkit/)
[![FastANI](https://img.shields.io/badge/FastANI-≥1.34.0-brightgreen.svg?style=flat)](https://github.com/ParBLiSS/FastANI)
[![KMA](https://img.shields.io/badge/KMA-≥1.4.0-brightgreen.svg?style=flat)](https://bitbucket.org/genomicepidemiology/kma/src/master/)

[![R](https://img.shields.io/badge/R-≥4.1.3-brightgreen.svg?style=flat)](https://www.r-project.org/)
[![python3](https://img.shields.io/badge/python-≥3.7-brightgreen.svg?style=flat)](https://www.python.org/downloads/)
