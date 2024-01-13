import argparse
import pandas as pd
import numpy as np
import re
import os

def calculate_flank_identity(input_file, outfile, metricsfile):
    # column data types
    dtype_dict = {
        'Sample': str,
        'Source': str,
        'Feature': str,
        'Start': int,
        'End': int,
        'Score': str, 
        'Strand': str,
        'Frame': str,
        'Attributes': str,
    }

    # load gff
    gff_df = pd.read_csv(input_file, sep='\t', header=None, names=dtype_dict.keys(), dtype=dtype_dict, low_memory=False)
    gene = re.sub(r'-flank_annotation\.gff$', '', os.path.basename(input_file))
    genome_dict = {}

    for index, row in gff_df.iterrows():
        genome_name = row['Sample']
        attributes = row['Attributes'].split(';')
        
        # Extract UniProtKB and ISfinder IDs
        ids = []
        for attr in attributes:
            if 'UniProtKB:' in attr:
                ids.append(attr.split('UniProtKB:')[1])
            elif 'ISfinder:' in attr:
                ids.append(attr.split('ISfinder:')[1])

        if genome_name in genome_dict:
            genome_dict[genome_name].extend(ids)
        else:
            genome_dict[genome_name] = ids

    genome_names = list(genome_dict.keys())
    num_samples = len(genome_names)

    #matrix full of NAs
    FI_matrix = np.full((num_samples, num_samples), np.nan)

    for i in range(num_samples):
        for j in range(i, num_samples):
            num1 = len(genome_dict[genome_names[i]])
            num2 = len(genome_dict[genome_names[j]])
            avg_num_genes = (num1 + num2) / 2

            if avg_num_genes < 2:
                continue

            intersection = len(set(genome_dict[genome_names[i]]) & set(genome_dict[genome_names[j]]))
            FI_score = (intersection / avg_num_genes) * 100
            FI_matrix[i, j] = FI_score
            FI_matrix[j, i] = FI_score

    # Create a df with sample names as row and column labels
    FI_df = pd.DataFrame(FI_matrix, index=genome_names, columns=genome_names)

    # save matrix
    FI_df.to_csv(outfile, sep='\t', index_label='Sample')
    
    # calculate metrics
    lower_triangle = FI_matrix[np.tril_indices(num_samples, k=-1)]
    median = np.nanmedian(lower_triangle)
    mean = np.nanmean(lower_triangle)

    with open(metricsfile, 'a') as metrics_file:
        metrics_file.write(f"{gene}\t{median}\t{mean}\n")

    print(f'{gene} complete')


parser = argparse.ArgumentParser(description='Calculate Flank identity matrix from GFF')
parser.add_argument('--gff', help='Path tp GFF file')
parser.add_argument('--outfile', help='Path to output file for Flank identity matrix')
parser.add_argument('--metricsfile', help='Path to output file for median and mean values')

args = parser.parse_args()
calculate_flank_identity(args.gff, args.outfile, args.metricsfile)
