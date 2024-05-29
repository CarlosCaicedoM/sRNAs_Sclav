#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 15:14:12 2024

@author: usuario
"""

import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt


#print(os.path.exists(gtf_guide))

# Define the folder containing the BAM files
bam_folder = '/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/results_rna_seq/S-clavuligerus_RNA-Seq-data_ALL/Bowtie2'

# Define the folder where you want to save the output GTF files
output_folder = '/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/results_rna_seq/S-clavuligerus_RNA-Seq-data_ALL'

# Path to the guide GTF file
gtf_guide = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/GCF_005519465.1_ASM551946v1_genomic.gff"



# Ensure the output folder exists
os.makedirs(output_folder, exist_ok=True)

# List all BAM files in the folder
bam_files = [f for f in os.listdir(bam_folder) if f.endswith('.bam')]

# List to store the generated GTF files
gtf_files = []

# Run StringTie for each BAM file
for bam_file in bam_files:
    bam_path = os.path.join(bam_folder, bam_file)
    output_path = os.path.join(output_folder, f"{os.path.splitext(bam_file)[0]}.gtf")
    gtf_files.append(output_path)
    
    # StringTie Command
    cmd = ['stringtie', bam_path, '-G', gtf_guide, 
           '-o', output_path, '-p', '20', '-m', '50', 
           '-s', '1', '-c', '1', '-g', '5']  
    """cmd = ['stringtie', bam_path, '-o', output_path,
           '-p', '20', '-m', '50', 
           '-s', '1', '-c', '1', '-g', '5']  """
    # Execute the command
    try:
        subprocess.run(cmd, check=True)
        print(f"StringTie procesado para {bam_file} exitosamente.")
    except subprocess.CalledProcessError as e:
        print(f"Error al procesar {bam_file}: {e}")

print("All BAM files have been processed.")


#%%  # Define the combined output GTF file
merged_gtf = os.path.join(output_folder, 'merged.gtf')

# Create a file listing all the GTF files
gtf_list_file = os.path.join(output_folder, 'gtf_list.txt')
with open(gtf_list_file, 'w') as f:
    for gtf_file in gtf_files:
        f.write(f"{gtf_file}\n")

# StringTie --merge command
#cmd_merge = ['stringtie', '--merge',  '-o', merged_gtf, gtf_list_file]
cmd_merge = ['stringtie', '--merge', '-G', gtf_guide, '-o', merged_gtf, gtf_list_file]
# Execute the merge command
try:
    subprocess.run(cmd_merge, check=True)
    print(f"GTF files merged successfully into {merged_gtf}.")
    print(f"GTF list file generated at {gtf_list_file}.")
except subprocess.CalledProcessError as e:
    print(f"Error merging GTF files: {e}")


# Output prefix for gffcompare
output_prefix = os.path.join(output_folder, 'gffcompare_output')

# gffcompare command with -i option
gffcompare_cmd = ['gffcompare', '-r', gtf_guide, '-o', output_prefix, '-i', gtf_list_file]

# Run gffcompare
try:
    subprocess.run(gffcompare_cmd, check=True)
    print("gffcompare completed successfully.")
except subprocess.CalledProcessError as e:
    print(f"Error running gffcompare: {e}")
    print(e.stderr.decode())




#%%  # Define the combined output GTF file without considering OMEGA results


def list_files_without_trimmed(directory):
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
    filtered_files = [f for f in files if "trimmed" not in f] 
    return filtered_files

files_non_size_selected = list_files_without_trimmed(output_folder)


def add_full_path(file_list, output_path):
    return [os.path.join(output_path, file) for file in file_list]

full_path_files_non_size_selected= add_full_path(files_non_size_selected, output_folder)






merged_gtf = os.path.join(output_folder, 'merged_non_size_selected.gtf')

# Create a file listing all the GTF files
gtf_list_file = os.path.join(output_folder, 'gtf_list.txt')
with open(gtf_list_file, 'w') as f:
    for gtf_file in full_path_files_non_size_selected:
        f.write(f"{gtf_file}\n")

# StringTie --merge command
#cmd_merge = ['stringtie', '--merge',  '-o', merged_gtf, gtf_list_file]
cmd_merge = ['stringtie', '--merge', '-G', gtf_guide, '-o', merged_gtf, gtf_list_file]
# Execute the merge command
try:
    subprocess.run(cmd_merge, check=True)
    print(f"GTF files merged successfully into {merged_gtf}.")
    print(f"GTF list file generated at {gtf_list_file}.")
except subprocess.CalledProcessError as e:
    print(f"Error merging GTF files: {e}")


# Output prefix for gffcompare
output_prefix = os.path.join(output_folder, 'gffcompare_output')

# gffcompare command with -i option
gffcompare_cmd = ['gffcompare', '-r', gtf_guide, '-o', output_prefix, '-i', gtf_list_file]

# Run gffcompare
try:
    subprocess.run(gffcompare_cmd, check=True)
    print("gffcompare completed successfully.")
except subprocess.CalledProcessError as e:
    print(f"Error running gffcompare: {e}")
    print(e.stderr.decode())






#%%  analyze merged assemblies

#open GFF file
output_folder = '/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/results_rna_seq/S-clavuligerus_RNA-Seq-data_ALL'

merged_gtf = os.path.join(output_folder, 
                          "Stringtie_GFFcompare_all_samples",
                          "merged.gtf")


gff_stringtie = pd.read_csv(merged_gtf, comment='#', 
                            sep = "\t", header=None)

gff_columns = ["seqname", "source", "feature", "start", "end", "score", 
                   "strand", "frame", "attribute" ]
gff_stringtie.columns = gff_columns

gff_stringtie_transcripts = gff_stringtie[gff_stringtie['feature'].str.startswith('transcript')]

StringTie_filtered = gff_stringtie_transcripts[~gff_stringtie_transcripts['attribute'].str.contains('gene-')]

#gff GFF compare

GFF_Compare_output = os.path.join(output_folder, 
                                  "Stringtie_GFFcompare_all_samples", 
                                  "gffcompare_output.combined.gtf")
GFFCompare = pd.read_csv(GFF_Compare_output, comment='#', 
                            sep = "\t", header=None)

gff_columns = ["seqname", "source", "feature", "start", "end", "score", 
                   "strand", "frame", "attribute" ]
GFFCompare.columns = gff_columns

GFFCompare_transcripts = GFFCompare[GFFCompare['feature'].str.startswith('transcript')]

GFFCompare_IGR = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "u"')]

GFFCompare_p = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "p"')]
GFFCompare_o = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "o"')]
GFFCompare_exact = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "="')]
GFFCompare_c = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "c"')]
GFFCompare_k = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "k"')]
GFFCompare_m = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "m"')]
GFFCompare_n = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "n"')]
GFFCompare_j = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "j"')]
GFFCompare_e = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "e"')]
GFFCompare_s = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "s"')]
GFFCompare_i = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "i"')]
GFFCompare_y = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "y"')]
GFFCompare_r = GFFCompare_transcripts[GFFCompare_transcripts['attribute'].str.contains('class_code "r"')]


def extract_num_samples(attribute):
    for item in attribute.split(';'):
        if 'num_samples' in item:
            return item.split('"')[1]

# Aplicar la función a la columna attributes
GFFCompare_IGR['num_samples'] = GFFCompare_IGR['attribute'].apply(extract_num_samples)

GFFCompare_p['num_samples'] = GFFCompare_p['attribute'].apply(extract_num_samples)

GFFCompare_o['num_samples'] = GFFCompare_o['attribute'].apply(extract_num_samples)

GFFCompare_exact['num_samples'] = GFFCompare_exact['attribute'].apply(extract_num_samples)
GFFCompare_exact['num_samples'] = pd.to_numeric(GFFCompare_exact['num_samples'])

GFFCompare_c['num_samples'] = GFFCompare_c['attribute'].apply(extract_num_samples)
GFFCompare_c['num_samples'] = pd.to_numeric(GFFCompare_c['num_samples'])



plt.hist(GFFCompare_exact['num_samples'], 
         bins=50, 
         edgecolor='black', 
         color = "darkorange")

# Añadir títulos y etiquetas
plt.xlabel('Number of samples')
plt.ylabel('Frequency')