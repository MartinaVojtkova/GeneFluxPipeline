# Snakefile
import csv
import numpy as np

configfile: "config.yaml"

#load samples from the TSV file
def load_samples(tsv_path):
    samples = {}
    with open(tsv_path) as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        next(reader)
        for row in reader:
            sample_name = row[0]  
            file_path = row[1]  
            samples[sample_name] = file_path
    return samples

# Load the TSV file specified in the config file
tsv_path = config["input_list"]
samples = load_samples(tsv_path)


rule all:
    input:
        expand("results/prokka_output/{genome}/{genome}.gff", genome=samples.keys()),
        "results/roary_results/clustered_proteins",
        "results/roary_results/gene_presence_absence.csv",
        "temp/gene_list.txt",
        "results/prokka_output/gene_db/gene_db.gff",
        "results/metrics/target_and_flank_median.dist",
        "results/metrics/flank_only_median.dist",
        "results/metrics/flank_identity.index",
        "results/metrics/flank_only_median.ani",
        "results/metrics/procrustes.results",
        "results/metrics/gene.metrics",
        "results/config.yaml"

rule dRep:
    input: 
       lambda wildcards: samples[wildcards.genome]
    output: 
        "results/dRep/data_tables/Mdb.csv" 
    params: 
        outdir = "results/dRep",
        pc= config["drep_pc"], 
        sc= config["drep_sc"]
    threads: 40
    conda: "environment.yaml"
    shell:
        """       
        dRep compare {params.outdir} -p {threads} -pa {params.pc} -sa {params.sc} -g {input}
        """

rule dRep_make_matrices: 
    input: 
        "results/dRep/data_tables/Mdb.csv"
    output: 
        dist = "results/dRep/matrices/genome_distance_matrix.txt" 
    params: 
        inlist = config["input_list"]
    conda: "environment.yaml"
    shell:
        """
        Rscript bin/dRep_make_matrices.R {input} {params.inlist} {output.dist} {output.sim}
        """

rule run_prokka:
    input:
        lambda wildcards: samples[wildcards.genome]
    output:
        "results/prokka_output/{genome}/{genome}.gff", 
        "results/prokka_output/{genome}/{genome}.fna"
    params: 
        genus =  config["genus"]
    threads: 40
    conda: "environment.yaml"
    shell:
        """
        prokka --cpus {threads} --kingdom Bacteria --genus {params.genus} --force --outdir results/prokka_output/{wildcards.genome} --prefix {wildcards.genome} --locustag {wildcards.genome} {input}      
        """

rule run_roary:
    input:
        expand("results/prokka_output/{genome}/{genome}.gff", genome=samples.keys())
    output:
        "results/roary_results/clustered_proteins",
        "results/roary_results/pan_genome_reference.fa", 
        "results/roary_results/gene_presence_absence.csv"
    params:
        outdir = "results/roary_results"
    threads: 40
    conda: "environment.yaml"
    shell:
        """
        if [ -d "{params.outdir}" ]; then
            rm -r "{params.outdir}"
        fi
        roary -p {threads} -f {params.outdir} -s -n -v -e -ap -g 70000 -o clustered_proteins {input}
        """

rule presence_absence_matrix:
    input:
        cluster_file = "results/roary_results/clustered_proteins",
        metadata = tsv_path
    output:
        "temp/presence_absence_matrix.tsv"
    params:
        outdir = "temp/",
    conda: "environment.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript bin/create_presence_absence_matrix.R {input.metadata} {input.cluster_file} {params.outdir} 
        """

rule make_gene_list:
    input:
        "temp/presence_absence_matrix.tsv"
    output:
        "temp/gene_list.txt"
    params:
        outdir = "temp/",
        slice_size= config["slice_size"]
    conda: "environment.yaml"
    shell:
        """
        Rscript bin/create_gene_list.R {input} {params.outdir} {params.slice_size}
        """

rule rename_reference_db:
   input:
      "results/roary_results/pan_genome_reference.fa"
   output:
      "temp/renamed_pan_genome_reference.fa"
   conda: "environment.yaml"
   shell:
    """
    awk '{{if($0 ~ /^>/) {{split($0,a," "); print ">"a[2]}} else print $0}}' {input} > {output}
    """

#extract gene_list.txt sequences 
rule create_flank_database:
   input:
     gene_list = "temp/gene_list.txt",
     reference = "temp/renamed_pan_genome_reference.fa"
   output:
     "results/flank_input/gene_db.fa"
   conda: "environment.yaml"
   shell:
     """
     seqkit grep -f {input.gene_list} -i {input.reference} > {output}
     """


rule annotate_gene_db: 
    input:
        "results/flank_input/gene_db.fa"
    output:
        "results/prokka_output/gene_db/gene_db.gff"
    conda: "environment.yaml"
    params:
       threads = 40
    shell:
        """
        prokka --cpus {params.threads} --kingdom Bacteria --genus Escherichia --force --outdir results/prokka_output/gene_db --prefix gene_db --locustag gene_db {input}      
        """

rule make_bedfiles:
    input:
        gff = expand("results/prokka_output/{genome}/{genome}.gff", genome=samples.keys()), 
        pa_file = "results/roary_results/gene_presence_absence.csv",
        gene_list = "temp/gene_list.txt"
    output:
        expand("results/bedfiles/{genome}/{genome}.bed", genome=samples.keys())
    params:
       threads = 40
    run:
        # Read gene list
        with open(input.gene_list, 'r') as file:
            gene_list = set(line.strip() for line in file)

        # Create gene dictionary
        gene_dict = {}
        with open(input.pa_file, 'r') as csvfile:
            reader = csv.reader(csvfile)
            headers = next(reader)
            first_genome_column = 14

            for row in reader:
                gene_name = row[0]
                if gene_name in gene_list:
                    for i in range(first_genome_column, len(headers)):
                        genome_name = headers[i]
                        gene_id = row[i]
                        if gene_id:
                            if genome_name not in gene_dict:
                                gene_dict[genome_name] = {}
                            gene_dict[genome_name][gene_id] = gene_name

        # Process each GFF file
        flank_distance = int(config["flank_length"])
        for gff_path in input.gff:
            filename_split = os.path.basename(gff_path).split('.')
            genome_name = '.'.join(filename_split[:2]) if len(filename_split) > 1 else filename_split[0]
            output_file_path = f"results/bedfiles/{genome_name}/{genome_name}.bed"

            with open(gff_path, 'r') as gff_file, open(output_file_path, 'w') as output_file:
                for line in gff_file:
                    if line.startswith('>'):
                        break
                    if line.startswith('#'):
                        continue

                    cols = line.strip().split('\t')
                    if len(cols) < 9:
                        continue

                    feat = cols[2]
                    if feat != 'CDS':
                        continue

                    attrib = cols[8]
                    attrib_dict = {item.split('=')[0]: item.split('=')[1] for item in attrib.split(';') if '=' in item}
                    gene_id = attrib_dict.get('ID')

                    if gene_id and gene_id in gene_dict.get(genome_name, {}):
                        start = int(cols[3]) - flank_distance
                        end = int(cols[4]) + flank_distance
                        if start < 0:
                            continue
                        gene_name = gene_dict[genome_name][gene_id]
                        output_line = f"{cols[0]}\t{start}\t{end}\t{gene_name}-{genome_name}\n"
                        output_file.write(output_line)
            print(f"Writing bedfile for {genome_name} complete")

             
rule bedtools_extract:
    input: 
        bedfile = "results/bedfiles/{genome}/{genome}.bed",
        fasta = "results/prokka_output/{genome}/{genome}.fna"
    output:
        fasta = "results/sequences/{genome}/{genome}_flanking_region.fasta"
    conda: "environment.yaml"
    shell:
         """
         bedtools getfasta -fi {input.fasta} -bed {input.bedfile} -name -fo {output}
         """

rule concatenate_sequences:
    input:
        expand("results/sequences/{genome}/{genome}_flanking_region.fasta", genome=samples.keys())
    output:
        "results/all_sequences_combined.fasta"
    shell:
        """
        rm results/sequences/E_coli_GCF_000005845.2/E_coli_GCF_000005845.2_flanking_region.fasta
        cat {input} > {output}
        """

checkpoint split_fasta_by_gene:
    input:
        "results/all_sequences_combined.fasta"
    output:
        directory("results/per_gene_sequences/target_and_flank")
    shell:
        """
        mkdir -p {output}
        awk '/^>/ {{
            split($$0, parts, "-");
            gene=substr(parts[1], 2);
            genome=parts[2];
            file="{output}/"gene".fasta";  
            print ">"genome > file; next;
        }} {{ print > file }}' {input}
        """

def flank_target_wildcards(wildcards):
    checkpoint_output = checkpoints.split_fasta_by_gene.get(**wildcards).output[0]
    return expand("results/per_gene_sequences/target_and_flank/{wildcard}.fasta", wildcard=glob_wildcards(os.path.join(checkpoint_output, "{wildcard}.fasta")).wildcard)

checkpoint mask_target_sequences:
    input:
        flank_target_wildcards
    output:
        directory("results/per_gene_sequences/flank_only")
    run:
        flank_length = int(config["flank_length"])
        output_dir = "results/per_gene_sequences/flank_only"

        # Create the output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Process each input file
        for input_file in input:
            gene = os.path.basename(input_file).split('.')[0]
            output_file = f"{output}/{gene}_flank_only.fasta"

            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                seq = ''
                header = None
                for line in infile:
                    if line.startswith('>'):
                        if header:
                            # Process the previous sequence
                            if len(seq) > 2 * flank_length:
                                NNN_seq = seq[:flank_length] + 'N' * (len(seq) - 2 * flank_length) + seq[-flank_length:]
                            else:
                                NNN_seq = seq  
                            outfile.write(header + NNN_seq + '\n')
                        header = line
                        seq = ''
                    else:
                        seq += line.strip()
                if header:
                    # Process the last sequence
                    if len(seq) > 2 * flank_length:
                        NNN_seq = seq[:flank_length] + 'N' * (len(seq) - 2 * flank_length) + seq[-flank_length:]
                    else:
                        NNN_seq = seq
                    outfile.write(header + NNN_seq + '\n')



rule index_flank_kma:
    input:
        target_and_flank="results/per_gene_sequences/target_and_flank/{gene}.fasta",
        flank_only="results/per_gene_sequences/flank_only/{gene}_flank_only.fasta"
    output:
        "results/kma_index/{gene}/target_and_flank.comp.b",
        "results/kma_index/{gene}/target_and_flank.length.b",
        "results/kma_index/{gene}/target_and_flank.name",
        "results/kma_index/{gene}/target_and_flank.seq.b",
        "results/kma_index/{gene}/flank_only.comp.b",
        "results/kma_index/{gene}/flank_only.length.b",
        "results/kma_index/{gene}/flank_only.name",
        "results/kma_index/{gene}/flank_only.seq.b"
    params:
        outfolder_target_and_flank=lambda wildcards: "results/kma_index/" +  wildcards.gene + "/target_and_flank",
        outfolder_flank_only=lambda wildcards: "results/kma_index/" +  wildcards.gene + "/flank_only",
        k=config["k-mer_size"]
    conda: "environment.yaml"
    shell:
        "kma index -k {params.k} -i {input.target_and_flank} -o {params.outfolder_target_and_flank};"
        "kma index -k {params.k} -i {input.flank_only} -o {params.outfolder_flank_only};"

rule flank_distance_matrix:
    input:
        ta="results/kma_index/{gene}/target_and_flank.comp.b",
        tb="results/kma_index/{gene}/target_and_flank.length.b",
        tc="results/kma_index/{gene}/target_and_flank.name",
        fa="results/kma_index/{gene}/flank_only.comp.b",
        fb="results/kma_index/{gene}/flank_only.length.b",
        fc="results/kma_index/{gene}/flank_only.name",
        t="results/kma_index/{gene}/target_and_flank.seq.b",
        f="results/kma_index/{gene}/flank_only.seq.b"
    output:
        t="results/kma_dist/{gene}/{gene}_target_and_flank.dist",
        f="results/kma_dist/{gene}/{gene}_flank_only.dist"
    params:
        t_db=lambda wildcards: "results/kma_index/" +  wildcards.gene + "/target_and_flank",
        f_db=lambda wildcards: "results/kma_index/" +  wildcards.gene + "/flank_only",
        dist_measure=1
     
    conda: "environment.yaml"
    shell:
        "kma dist -t_db {params.t_db} -o {output.t} -d {params.dist_measure};"
        "kma dist -t_db {params.f_db} -o {output.f} -d {params.dist_measure};"

def kma_tf(wildcards):
    checkpoint_output = checkpoints.mask_target_sequences.get(**wildcards).output[0]
    return expand("results/kma_dist/{gene}/{gene}_target_and_flank.dist", gene=glob_wildcards(os.path.join(checkpoint_output, "{gene}_flank_only.fasta")).gene)

def kma_f(wildcards):
    checkpoint_output = checkpoints.mask_target_sequences.get(**wildcards).output[0]
    return expand("results/kma_dist/{gene}/{gene}_flank_only.dist", gene=glob_wildcards(os.path.join(checkpoint_output, "{gene}_flank_only.fasta")).gene)


rule median_target_and_flank_dist:
    input:
        kma_tf
    output:
        t_dist = "results/metrics/target_and_flank_median.dist"
    run:
        flank_length = int(config["flank_length"])
        output_dir = "results/metrics"
        os.makedirs(output_dir, exist_ok=True)

        with open(output.t_dist, 'w') as output_file:  
            
            output_file.write("Gene\tDistance_median\n")

            for filename in input:
                gene = re.sub(r'_target_and_flank\.dist$', '', os.path.basename(filename))

                with open(filename, 'r') as file:
                    lines = file.readlines()

                size = int(lines[0].strip())
                matrix = np.zeros((size, size))

                for i, line in enumerate(lines[1:]):
                    values = line.split()[1:]  # Skip the first element (row name)
                    numeric_values = [int(val) if val.isdigit() else 2 * flank_length for val in values]
                    matrix[i, :len(numeric_values)] = numeric_values

                for i in range(size):
                    for j in range(i):
                        matrix[j, i] = matrix[i, j]

              
                lower_triangle = matrix[np.tril_indices(size, k=-1)]
                median = np.median(lower_triangle)

                print(f"Added: {gene}")
                output_file.write(f"{gene}\t{median}\n")

rule median_flank_only_dist:
    input:
        kma_f
    output:
        t_dist = "results/metrics/flank_only_median.dist"
    run:
        flank_length = int(config["flank_length"])
        output_dir = "results/metrics"
        os.makedirs(output_dir, exist_ok=True)

        with open(output.t_dist, 'w') as output_file:  
            
            output_file.write("Gene\tDistance_median\n")

            for filename in input:
                gene = re.sub(r'_flank_only\.dist$', '', os.path.basename(filename))
                print(f"processing {gene}")
                with open(filename, 'r') as file:
                    lines = file.readlines()

                size = int(lines[0].strip())
                matrix = np.zeros((size, size))

                for i, line in enumerate(lines[1:]):
                    values = line.split()[1:]  
                    numeric_values = [int(val) if val.isdigit() else 2 * flank_length for val in values]
                    matrix[i, :len(numeric_values)] = numeric_values

                for i in range(size):
                    for j in range(i):
                        matrix[j, i] = matrix[i, j]

              
                lower_triangle = matrix[np.tril_indices(size, k=-1)]
                median = np.median(lower_triangle)

                print(f"Added: {gene}")
                output_file.write(f"{gene}\t{median}\n")


rule get_flank_annotations_per_genome:
    input: 
        bed = "results/bedfiles/{genome}/{genome}.bed",
        gff = "results/prokka_output/{genome}/{genome}.gff"
    output: 
        directory("results/flank_annotations/{genome}")
    conda: "environment.yaml"
    shell:
        """
        mkdir -p {output}
        python3 bin/bed_to_flank_annotation.py --gff {input.gff} --bed {input.bed} --outdir {output}
        """

 
checkpoint merge_gene_annotations:
    input:
        genome_dirs = expand("results/flank_annotations/{genome}", genome=samples.keys())
    output:
        merged_dir = directory("results/merged_annotations")
    conda:
        "environment.yaml"
    shell:
       """
        mkdir -p {output.merged_dir}
        for genome_dir in {input.genome_dirs}; do
            for gff_file in $genome_dir/*.gff; do
                gene=$(basename $gff_file .gff)
                output_file="{output.merged_dir}/${{gene}}-flank_annotation.gff"
                if [ ! -f $output_file ]; then
                    cat $gff_file > $output_file
                else
                    cat $gff_file >> $output_file
                fi
            done
        done
        """

def merged_annot_files(wildcards):
    checkpoint_output = checkpoints.merge_gene_annotations.get(**wildcards).output[0]
    return expand("results/merged_annotations/{gene}-flank_annotation.gff", gene=glob_wildcards(os.path.join(checkpoint_output, "{gene}_flank_annotation.gff")).gene)


rule separate_gene_fasta:
    input:
        "results/per_gene_sequences/flank_only/{gene}_flank_only.fasta"
    output:
        directory("results/fastani_input/sequences/{gene}"),
    conda:
        "environment.yaml"
    shell:
        """
        mkdir -p {output}
        awk '/^>/ {{if (outfile!=""){{close(outfile)}}; outfile="{output}/" substr($0,2) ".fasta"}} {{print >> outfile}}' {input}
        """



def list_input(wildcards):
    checkpoint_output = checkpoints.mask_target_sequences.get(**wildcards).output[0]
    return expand("results/fastani_input/sequences/{gene}", gene=glob_wildcards(os.path.join(checkpoint_output, "{gene}_flank_only.fasta")).gene)

checkpoint create_lists:
    input:
        list_input
    output:
        directory("results/fastani_input/lists")
    conda: "environment.yaml"
    shell:
        """
        mkdir -p {output}
        for folder in {input}; do
            gene=$(basename "$folder")
            find "$folder" -type f -name '*.fasta' | sort > "{output}/${{gene}}_fasta_list.txt"
        done
        """

rule fastani:
    input: 
        "results/fastani_input/lists/{gene}_fasta_list.txt"
    output:
        outfile = "results/fastani_results/{gene}/{gene}_ANI", 
        matrix = "results/fastani_results/{gene}/{gene}_ANI.matrix"
    conda: "environment.yaml"
    shell:
        """
        fastANI --ql {input} --rl {input} -o {output.outfile} --matrix
        """

def ANI_files(wildcards):
    checkpoint_output = checkpoints.create_lists.get(**wildcards).output[0]
    return expand("results/fastani_results/{gene}/{gene}_ANI", gene=glob_wildcards(os.path.join(checkpoint_output, "{gene}_fasta_list.txt")).gene)

def ANI_matrices(wildcards):
    checkpoint_output = checkpoints.create_lists.get(**wildcards).output[0]
    return expand("results/fastani_results/{gene}/{gene}_ANI.matrix", gene=glob_wildcards(os.path.join(checkpoint_output, "{gene}_fasta_list.txt")).gene)

rule flank_index_calculation:
    input: 
        gff = "results/merged_annotations/{gene}-flank_annotation.gff"
    output: 
        matrix = "results/flank_identity/{gene}_flank_identity.matrix",
        metrics = "temp/metrics/{gene}_FI_metrics.txt"
    params:
        matrixdir = "results/flank_identity",
        metricsdir = "temp/metrics"
    conda: "environment.yaml"
    shell:
        """
        mkdir -p {params.matrixdir}
        mkdir -p {params.metricsdir}
        python3 bin/flank_identity.py --gff {input.gff} --outfile {output.matrix} --metricsfile {output.metrics}
        """

def matrix_input(wildcards):
    checkpoint_output = checkpoints.mask_target_sequences.get(**wildcards).output[0]
    return expand("results/flank_identity/{gene}_flank_identity.matrix", gene=glob_wildcards(os.path.join(checkpoint_output, "{gene}_flank_only.fasta")).gene)

def metrics_input(wildcards):
    checkpoint_output = checkpoints.mask_target_sequences.get(**wildcards).output[0]
    return expand("temp/metrics/{gene}_FI_metrics.txt", gene=glob_wildcards(os.path.join(checkpoint_output, "{gene}_flank_only.fasta")).gene)


rule merge_FI_metrics:
    input: 
        matrix_input,
        metrics_input
    output: 
        "results/metrics/flank_identity.index"
    shell:
        """
        echo -e "Gene\tFlank_index_median\tFlank_index_mean" > {output} && cat temp/metrics/*_FI_metrics.txt >> {output}
        """

rule median_flank_only_ani:
    input:
        expand("results/fastani_results/{gene}/{gene}_ANI.matrix", gene=glob_wildcards("results/fastani_input/lists/{gene}_fasta_list.txt").gene)
    output:
        t_dist = "results/metrics/flank_only_median.ani"
    run:
        output_dir = "results/metrics"
        os.makedirs(output_dir, exist_ok=True)

        with open(output.t_dist, 'w') as output_file:  
            
            output_file.write("Gene\tANI_median\n")

            for filename in input:
                gene = re.sub(r'_ANI\.matrix$', '', os.path.basename(filename))

                with open(filename, 'r') as file:
                    lines = file.readlines()

                size = int(lines[0].strip())
                matrix = np.zeros((size, size))
                

                for i, line in enumerate(lines[1:]):
                    values = line.split()[1:]  # Skip the first element (row name)
                    numeric_values = [float(val) if val != 'NA' else np.nan for val in values]
                    matrix[i, :len(numeric_values)] = numeric_values

                for i in range(size):
                    for j in range(i):
                        matrix[j, i] = matrix[i, j]

                lower_triangle = matrix[np.tril_indices(size, k=-1)]
                median = np.nanmedian(lower_triangle)

                print(f"Added: {gene}")
                output_file.write(f"{gene}\t{median}\n")


rule procrustes:
    input:
        kma = "results/kma_dist/{gene}/{gene}_flank_only.dist",
        dRep = "results/dRep/matrices/genome_distance_matrix.txt"
    output:
        "results/procrustes/temp/{gene}_proc.metrics"
    params: 
        outdir = "results/procrustes/temp"
    conda: "environment.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript bin/procrustes.R {input.dRep} {input.kma} {output}         
        """

def proc_column_input(wildcards):
    checkpoint_output = checkpoints.mask_target_sequences.get(**wildcards).output[0]
    return expand("results/procrustes/temp/{gene}_proc.metrics", gene=glob_wildcards(os.path.join(checkpoint_output, "{gene}_flank_only.fasta")).gene)

rule merge_proc_results:
    input: 
        proc_column_input
    output: 
        "results/metrics/procrustes.results"
    shell:
        """
        echo -e "Gene\tProc_Correlation\tProc_Signif" > {output} && cat results/procrustes/temp/*_proc.metrics >> {output}
        """


rule mantel:
    input:
        ani = "results/fastani_results/{gene}/{gene}_ANI.matrix",
        dRep = "results/dRep/matrices/genome_similarity_matrix.txt"
    output:
        "results/mantel/temp/{gene}_mantel.metrics"
    params: 
        outdir = "results/mantel/temp"
    conda: "environment.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript bin/mantel.R {input.dRep} {input.ani} {output}         
        """

def mantel_column_input(wildcards):
    checkpoint_output = checkpoints.mask_target_sequences.get(**wildcards).output[0]
    return expand("results/mantel/temp/{gene}_mantel.metrics", gene=glob_wildcards(os.path.join(checkpoint_output, "{gene}_flank_only.fasta")).gene)


rule merge_man_results:
    input: 
        mantel_column_input
    output: 
        "results/metrics/mantel.results"
    shell:
        """
        echo -e "Gene\tMan_Correlation\tMan_Signif" > {output} && cat results/procrustes/temp/*_proc.metrics >> {output}
        """


rule merge_metrics:
    input:
        FI="results/metrics/flank_identity.index",
        KMA="results/metrics/flank_only_median.dist",
        ANI="results/metrics/flank_only_median.ani", 
        PROC="results/metrics/procrustes.results", 
        MAN = "results/metrics/mantel.results"
    output:
        "results/metrics/gene.metrics"
    run:
        import pandas as pd

        fi = pd.read_csv(input.FI, sep='\t')
        kma = pd.read_csv(input.KMA, sep='\t')
        ani = pd.read_csv(input.ANI, sep='\t')
        proc = pd.read_csv(input.PROC, sep='\t')
        man=pd.read_csv(input.MAN, sep='\t')
        merged_df = fi.merge(kma, on='Gene', how='outer').merge(ani, on='Gene', how='outer').merge(proc, on='Gene', how='outer').merge(man, on='Gene', how='outer')
        merged_df.fillna('nan', inplace=True)
        merged_df.to_csv(str(output[0]), sep='\t', index=False)


rule end:
    input: 
        merged_annot_files, 
        ANI_files
    output:
        config = "results/config.yaml"
    conda: "environment.yaml"
    shell: 
        """
        cp config.yaml {output.config}
        """

onsuccess:
    print("This is the end of the pipeline")