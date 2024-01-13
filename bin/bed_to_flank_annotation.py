import argparse

def parse_bed(bed_file):
    bed_lines = []
    with open(bed_file, 'r') as file:
        for line in file:
            cols = line.strip().split('\t')
            if len(cols) < 4:
                continue
            genome, start, end, info = cols
            gene_name = info.split('-')[0]
            bed_lines.append((genome, int(start), int(end), gene_name))
    return bed_lines

def database_gff(gff_file):
    gff_dict = {}
    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith('#') or line.startswith('>') or not line.strip():
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            genome = cols[0]
            start, end = int(cols[3]), int(cols[4])
            if genome not in gff_dict:
                gff_dict[genome] = []
            gff_dict[genome].append((start, end, line))
    return gff_dict

def filter_gff(bed_single, gff_dict):
    genome, start, end, gene_name = bed_single
    if genome not in gff_dict:
        return []
    
    return [line for gff_start, gff_end, line in gff_dict[genome] 
            if gff_start <= end and gff_end >= start]

def write_single_gffs(flank_annot, gene_name, outdir):
    output_file_path = f"{outdir}/{gene_name}.gff"
    with open(output_file_path, 'w') as file:
        for gene in flank_annot:
            file.write(gene)


def main(gff_file, bed_file, outdir):
    bed_lines = parse_bed(bed_file)
    gff_dict = database_gff(gff_file)

    for bed_single in bed_lines:
        genome, start, end, gene_name = bed_single
        flank_annot = filter_gff(bed_single, gff_dict)
        write_single_gffs(flank_annot, gene_name, outdir)


input = argparse.ArgumentParser(description='BED to flank GFFs')
input.add_argument('--gff', type=str, help='Path to GFF file')
input.add_argument('--bed', type=str, help='Path to BED file')
input.add_argument('--outdir', type=str, help='Path to output folder for gene GFFs')

args = input.parse_args()

main(args.gff, args.bed, args.outdir)

print(f"GFF flank files for {args.gff} generated")
