#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 09:44:08 2023

@author: usuario
"""



#conda install -c bioconda fastqc
#sudo apt install sra-toolkit
#conda install -c bioconda -c conda-forge multiqc
#sudo apt install cutadapt
#conda install -c bioconda fastp
#sudo apt install bwa
#sudo apt install bowtie2
#pip install pysam
#sudo apt install samtools
#pip install HTSeq
#pip install pydeseq2
#conda install -c conda-forge scanpy python-igraph leidenalg
#pip install gseapy
#pip install goatools
#pip install pybedtools
# conda install -c bioconda bbmap




#%% Import libraries required

from IPython import get_ipython
get_ipython().magic('reset -sf')


import subprocess
import os
import shutil
import pandas as pd
import re
import pickle as pkl
import pysam
import pybedtools
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import scanpy as sc
import gseapy as gp
from gseapy.plot import gseaplot
import numpy as np
import matplotlib.pyplot as plt
import concurrent.futures
import multiprocessing


#%% Input data


my_experiment = "S-clav_RNA-seq_sRNas_OMEGA"
reference_genome_fna = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/GCF_005519465.1_ASM551946v1_genomic.fna"
reference_genome_gtf = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/GCF_005519465.1_ASM551946v1_genomic.gtf"
reference_genome = "SCLAV"

reference_chr_fasta = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/NZ_CP027858.1/NZ_CP027858.1.fasta"
reference_plasmid_fasta = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/NZ_CP027859.1/NZ_CP027859.1.fasta"
reference_chr = "SCLAV_CHR"
reference_plasmid = "SCLAV_PLASMID"


cur_dir = os.getcwd()
#path_rawdata = os.path.join(cur_dir, "raw_data")
path_rawdata = os.path.join(cur_dir, "raw_data", "RNA-seq_OMEGA")


fastq_files = [file for file in os.listdir(path_rawdata) if file.endswith(".gz") ]
fastq_files_path = []
for i in fastq_files:
    fastq_files_path.append(os.path.join(path_rawdata, i))



#%% Initial quality control of the data with fastqc
my_experiment = "S-clav_RNA-seq_sRNas_OMEGA"

try:
    os.mkdir("results_rna_seq")
except OSError as error:
    print(error)  

try:
    os.mkdir(os.path.join("results_rna_seq", my_experiment))
except OSError as error:
    print(error)

    
    
    

    
fastqc_outdir=os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastqc")
try:
    os.mkdir(fastqc_outdir)
except OSError as error:
    print(error) 

#As Java version of PC is nos suitable for FASTQC I decided to create a conda 
#environment to run fastqc
def run_in_conda_env(conda_env_name, command_list):
    command_str = ' '.join(command_list)
    # Use conda run to execute the command in the specified environment
    full_command = ["conda", "run", "-n", conda_env_name] + command_list
    subprocess.run(full_command)

for i in fastq_files:
    fastQC_command = ["fastqc", os.path.join(path_rawdata, i), "-o", fastqc_outdir, "-t", "20"]
    run_in_conda_env("FASTQC", fastQC_command)

#Summarize all the fastqc results with multiQC
def MultiQC(inputdir, outdir):
    try:
        os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, outdir))
    except OSError as error:
        print(error)
        
    multiQC_command = ["multiqc", inputdir]
    multiQC_command.append("-n")
    multiQC_command.append(my_experiment)
    multiQC_command.append("-o")
    multiQC_command.append(os.path.join(cur_dir, "results_rna_seq", my_experiment,  outdir))
    multiQC_command.append("-f")
    multiQC_command.append("--profile-runtime") 
    subprocess.call(multiQC_command)

MultiQC(fastqc_outdir, "multiqc")

#%% Quality control and filtering by fastp

fastp_dir =os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp")
try:
  os.mkdir(fastp_dir)
except OSError as error:
    print(error)

#if single-end

for i,j  in enumerate(fastq_files):
    fastp_command=["fastp", 
                   "-i", fastq_files_path[i], "-o",  
                   os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+j),
                   "-j", j+".json", "-h", j+".html", "-w", "16", 
                   "--cut_right", "--trim_poly_g", "--trim_poly_x"]

    subprocess.call(fastp_command)

#If paired-end
experiments = []
for file in fastq_files:
    parts = file.split('_R')
    name = parts[0]
    experiments.append(name)
    
experiments = set(experiments) 

for i,j  in enumerate(experiments):
    fastp_command=["fastp", 
                   "-i", path_rawdata+"/"+j+"_R1_001.fastq.gz",
                   "-I",  path_rawdata+"/"+j+"_R2_001.fastq.gz",          
                   "-o", os.path.join(cur_dir, "results_rna_seq",
                                my_experiment, "fastp", "trimmed_"+j+"_R1.fastq.gz"), 
                   "-O", os.path.join(cur_dir, "results_rna_seq",
                                my_experiment, "fastp", "trimmed_"+j+"_R2.fastq.gz"),
                   "-j", j+".json", "-h", j+".html", "-w", "16", 
                   "--cut_right", "--trim_poly_g"]

    subprocess.call(fastp_command)

#Move the reports to the fastp folder
reports = [rep for rep in os.listdir(os.curdir) if (rep.endswith(".html") or rep.endswith(".json"))]  
for file_name in reports:
    shutil.move(file_name, os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp"))





#%% Fastqc after trimming

fastq_files_trimmed = [file for file in os.listdir(fastp_dir) if file.endswith(".gz") ]
fastQC_trimmed=os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastQC_trimmed")
try:
    os.mkdir(fastQC_trimmed)
except OSError as error:
    print(error) 

#As Java version of PC is nos suitable for FASTQC I decided to create a conda 
#environment to run fastqc


for i in fastq_files_trimmed:
    fastQC_command_fastp = ["fastqc", os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", i), "-o", fastQC_trimmed, "-t", "20"]
    run_in_conda_env("FASTQC", fastQC_command_fastp)


#Summarize all the fastqc results with multiQC
MultiQC(fastQC_trimmed, "multiqc_fastp")


#%%  Remove ribosomal RNA sequences


fastq_files_trimmed = [file for file in os.listdir(fastp_dir) if file.endswith(".gz") ]
fastq_rRNA_removed=os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed")

try:
    os.mkdir(fastq_rRNA_removed)
except OSError as error:
    print(error) 


for i, j  in enumerate(experiments):
    try:
        os.mkdir(j)
    except OSError as error:
        print(error) 
 
    sortmeRNA_command = ["sortmerna",
                         "--ref", "rRNA_databases/silva-bac-16s-id90.fasta", 
                         "--ref", "rRNA_databases/silva-bac-23s-id98.fasta", 
                         "--ref", "rRNA_databases/rfam-5s-database-id98.fasta", 
                         "--ref", "rRNA_databases/rfam-5.8s-database-id98.fasta", 
                         "--reads", os.path.join(cur_dir, "results_rna_seq",
                                      my_experiment, "fastp", "trimmed_"+j+"_R1.fastq.gz"), 
                         "--reads", os.path.join(cur_dir, "results_rna_seq",
                                      my_experiment, "fastp", "trimmed_"+j+"_R2.fastq.gz"),
                         "--other",  os.path.join(fastq_rRNA_removed, j), 
                         "--fastx", "--threads",  "20", "--paired_in", 
                         "--workdir", j] 

    subprocess.call(sortmeRNA_command)
    
#%% reformat the output of sortmeRNA to have again the forward and reverse fastq files
#as indepedent files

fastq_files_rRNA_removed= [file for file in os.listdir(fastq_rRNA_removed) if file.endswith(".gz") ]


for i, j  in enumerate(fastq_files_rRNA_removed):
   reformat_command = ["reformat.sh", 
                        "in="+os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed", j), 
                        "out1="+os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed", j.split(".fq")[0]+"_1.fastq.gz"),
                        "out2="+os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed", j.split(".fq")[0]+"_2.fastq.gz")]
   subprocess.call(reformat_command)


#%% Removeinterleaved files

character = '.fastq.'
for filename in os.listdir(fastq_rRNA_removed):
    if character not in filename:
        file_path = os.path.join(fastq_rRNA_removed, filename)
        # Remover el archivo
        os.remove(file_path)


#%% Fastqc after rRNA removal

fastq_files_rRNA_removed= [file for file in os.listdir(fastq_rRNA_removed) if file.endswith(".gz") ]
fastQC_rRNA_removed=os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastQC_rRNA_removed")
try:
    os.mkdir(fastQC_rRNA_removed)
except OSError as error:
    print(error) 

#As Java version of PC is nos suitable for FASTQC I decided to create a conda 
#environment to run fastqc

for i in fastq_files_rRNA_removed:
    fastQC_command_rRNA_removed = ["fastqc", os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed", i), "-o", fastQC_rRNA_removed, "-t", "20"]
    run_in_conda_env("FASTQC", fastQC_command_rRNA_removed)

#Summarize all the fastqc results with multiQC
MultiQC(fastQC_rRNA_removed, "multiqc_fastp_rRNA_removal")



#%% Align sequences to the reference genome  Bowtie2

# generate the index genome file
# Comando Bowtie2 para crear la base de datos de referencia
bowtie2_build_cmd = ["bowtie2-build", reference_genome_fna, reference_genome]
subprocess.run(bowtie2_build_cmd, check=True)
fastq_filtered_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp")
fastq_filtered = [fastq for fastq in os.listdir(fastq_filtered_folder) if fastq.endswith(".gz") ]  


#SE data
try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2"))
except OSError as error:
    print(error)
    
for i in fastq_filtered:  
    fastq_file = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", i)
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "Bowtie2", i)  +".sam"
    bowtie2_align_cmd = ["bowtie2", "--local", "-p", "18", "-x",  
                         reference_genome, "-U", fastq_file,  "-S",  output_sam]
    subprocess.run(bowtie2_align_cmd, check=True)


# Alignment for PE data WITHOUT rRNA removal 
try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_fastp_local"))
except OSError as error:
    print(error)

for i in experiments:  
    # Comamand for 
    input_fastq_r1 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_R1.fastq.gz")
    input_fastq_r2 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_R2.fastq.gz")
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "Bowtie2_fastp_local", i)  +".sam"
    bowtie2_align_cmd = ["bowtie2", "--local", "-p", "18", "-x", reference_genome, "-1", input_fastq_r1, 
                         "-2", input_fastq_r2, "-S", output_sam]
    subprocess.run(bowtie2_align_cmd, check=True)


# Alignment for PE data WITH rRNA removal 
try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_rRNA_removed"))
except OSError as error:
    print(error)

for i in experiments:  
    # Comamand for 
    input_fastq_r1 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed", i+"_1.fastq.gz")
    input_fastq_r2 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed", i+"_2.fastq.gz")
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "Bowtie2_rRNA_removed", i)  +".sam"
    bowtie2_align_cmd = ["bowtie2", "-p", "18", "-x", reference_genome, "-1", input_fastq_r1, 
                         "-2", input_fastq_r2, "-S", output_sam]
    subprocess.run(bowtie2_align_cmd, check=True)


#%% Convert SAM files to BAM files 

def convert_sam_to_bam(sam_file, bam_file):
    # Open SAM file in reading mode
    sam = pysam.AlignmentFile(sam_file, "r")

    # Open BAM file in writing mode
    #"When using 'wbu' as the opening mode for the BAM file, 
    #it will be written to the file uncompressed."
    bam = pysam.AlignmentFile(bam_file, "wb", header=sam.header)

    # Read each register in SAM file and write in BAM file
    for read in sam:
        bam.write(read)

    # Close the files
    sam.close()
    bam.close()
    
# Alignment for PE data WITHOUT rRNA removal 
alignment_sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_fastp_local")
alignment_sam = [sam for sam in os.listdir(alignment_sam_folder) if sam.endswith(".sam") ]  
for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    output_bam = os.path.join(alignment_sam_folder, os.path.splitext(alignment)[0]+".bam")
    convert_sam_to_bam(input_sam, output_bam)

# Alignment for PE data WITH rRNA removal 
alignment_sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_rRNA_removed")
alignment_sam = [sam for sam in os.listdir(alignment_sam_folder) if sam.endswith(".sam") ]  
for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    output_bam = os.path.join(alignment_sam_folder, os.path.splitext(alignment)[0]+".bam")
    convert_sam_to_bam(input_sam, output_bam)


#%% Convert unsorted to sorted BAM files

def check_bam_order(file_path):
    try:
        with pysam.AlignmentFile(file_path, 'rb') as bamfile:
            # Get the BAM file header
            header = bamfile.header
            # Check if the BAM file is sorted
            is_sorted = 'HD' in header and 'SO' in header['HD']
            return is_sorted
    except Exception as e:
        print(f"Error opening the BAM file: {e}")
        return False
    
def sorted_bam(alignment_sam_folder):    
    alignment_bam = [bam for bam in os.listdir(alignment_sam_folder) if bam.endswith(".bam") ]  
    for bam in alignment_bam: 
        input_bam = os.path.join(alignment_sam_folder, bam)    
        # Output file CRAM
        output_bam = os.path.join(alignment_sam_folder, "sorted_"+bam)
        
        # Command 1: samtools fixmate
        fixmate_cmd = ["samtools", "fixmate", "-m", input_bam, "-"]
        
        # Command 2: samtools sort
        sort_cmd = ["samtools", "sort", "-@20", "-T", os.path.join(alignment_sam_folder, "sort_tmp"), "-"]
        
        
        # Comando 3: samtools markdup
        markdup_cmd = ["samtools", "markdup", "-O", "bam", "-@20", "-", output_bam]
    
       
        # Ejecutar los comandos en secuencia
        try:
            # Ejecutar comando 1
            fixmate_process = subprocess.Popen(fixmate_cmd, stdout=subprocess.PIPE)
            
            # Ejecutar comando 2 y establecer la entrada del proceso anterior como la salida
            sort_process = subprocess.Popen(sort_cmd, stdin=fixmate_process.stdout, stdout=subprocess.PIPE)
            fixmate_process.stdout.close()
            
            # Ejecutar comando 3 y establecer la entrada del proceso anterior como la salida
            markdup_process = subprocess.Popen(markdup_cmd, stdin=sort_process.stdout)
            sort_process.stdout.close()
            
            # Esperar a que se completen todos los procesos
            fixmate_process.wait()
            sort_process.wait()
            markdup_process.wait()
        
        except subprocess.CalledProcessError as e:
            print("Error:", e)
sorted_bam(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_fastp_local"))
sorted_bam(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_rRNA_removed"))

#%% Check if the rRNA was succesfully removed by counting the number of 
#reads aligned to rRNA genes

reference_genome_gtf = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/GCF_005519465.1_ASM551946v1_genomic.gtf"
reference_genome_bed = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/GCF_005519465.1_ASM551946v1_genomic.bed"
alignment_bam_sorted = [bam for bam in os.listdir(alignment_sam_folder) if bam.startswith("sorted")]  

def convert_gff_to_bed(gff_file, bed_file):
    cmd = f"gtf2bed < {gff_file} > {bed_file}"
    subprocess.run(cmd, shell=True)
# Usar la función para convertir un archivo GFF a BED
convert_gff_to_bed(reference_genome_gtf, reference_genome_bed)


reference_genome_bed_rRNA = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/GCF_005519465.1_ASM551946v1_rRNA.bed"
command = f"grep 'gbkey \"rRNA\"' {reference_genome_bed} > {reference_genome_bed_rRNA}"
subprocess.run(command, shell=True)

#access a BAM index file
def bam_index_file(alignment_sam_folder, alignment_bam_sorted):
    for i in alignment_bam_sorted:
        bam_file = os.path.join(alignment_sam_folder, i)
        command_index = ["samtools", "index", bam_file]    
        subprocess.run(command_index)


for i in alignment_bam_sorted:

    bam_file = pysam.AlignmentFile(os.path.join(alignment_sam_folder, i), 'rb')
    bed_file = pybedtools.BedTool(reference_genome_bed_rRNA)
    
    # Filtrar lecturas que se alinean a regiones de rRNA
    filtered_reads = bed_file.intersect(os.path.join(alignment_sam_folder, i), u=True)
    
    # Contar lecturas en archivos BAM
    total_reads = sum(1 for read in bam_file)
    rRNA_reads = sum(1 for read in filtered_reads)
    
    # Calcular porcentaje
    percentage = 100.0 * rRNA_reads / total_reads
    
    print(f'El porcentaje de lecturas que se alinearon con rRNA es: {percentage}%')
    
    # Cerrar el archivo BAM
    bam_file.close()


#%%  Delete sam files


def delete_files_aln(folder, suffix):
    files = os.listdir(folder)
    for file in files:
        path_file = os.path.join(folder, file)
        if file.endswith(suffix) or (file.endswith(".bam") and not file.startswith("sorted")):
            os.remove(path_file)


delete_files_aln(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_fastp_local"), ".sam")
delete_files_aln(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_rRNA_removed"), ".sam")

#%% Align sequences to the reference genome with BWA

try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "BWA"))
except OSError as error:
    print(error)


# generate the index genome file
command = ['bwa', 'index', reference_genome_fna]
subprocess.run(command, check=True)


fastq_filtered_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp")
fastq_filtered = [fastq for fastq in os.listdir(fastq_filtered_folder) if fastq.endswith(".gz") ]  

#SE data
for i in fastq_filtered:  
    # Comamand for BWA
    fastq_file = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", i)
    command = ["bwa", "mem", reference_genome_fna, fastq_file, "-t", "18" ]
    # redirect the output to the SAM file
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "BWA", i)  +".sam"
    
    with open(output_sam, "w") as output_file:
        subprocess.run(command, stdout=output_file, check=True)
    


#PE data

for i in experiments:  
    # Comamand for BWA
    fastq_file_1 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_R1.fastq.gz")
    fastq_file_2 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_R2.fastq.gz")
    command_BWA = ["bwa", "mem", "-t", "18", reference_genome_fna, fastq_file_1, fastq_file_2]
    # redirect the output to the SAM file
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "BWA", i)  +".sam"
    
    with open(output_sam, "w") as output_file:
        subprocess.run(command_BWA, stdout=output_file, check=True)


#Convert SAM files to BAM files 
alignment_sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "BWA")
alignment_sam = [sam for sam in os.listdir(alignment_sam_folder) if sam.endswith(".sam") ]  
for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    output_bam = os.path.join(alignment_sam_folder, os.path.splitext(alignment)[0]+".bam")
    convert_sam_to_bam(input_sam, output_bam)




#Convert unsorted to sorted BAM files

for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    is_sorted = check_bam_order(input_sam)

if is_sorted:
    print("The BAM file is already sorted.")
else:
    sorted_bam(os.path.join(cur_dir, "results_rna_seq", my_experiment, "BWA"))

#Delete sam files and unsorted bam
delete_files_aln(os.path.join(cur_dir, "results_rna_seq", my_experiment, "BWA"), ".sam")

#%% Align sequences to the reference genome  Bowtie 

# generate the index genome file
# Comando Bowtie para crear la base de datos de referencia
bowtie_build_cmd = ["bowtie-build", reference_genome_fna, reference_genome]
subprocess.run(bowtie_build_cmd, check=True)

fastq_filtered_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp")
fastq_filtered = [fastq for fastq in os.listdir(fastq_filtered_folder) if fastq.endswith(".gz") ]  


# Alignment for PE data WITHOUT rRNA removal 
try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie"))
except OSError as error:
    print(error)

for i in experiments:  
    # Comamand for 
    input_fastq_r1 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_R1.fastq.gz")
    input_fastq_r2 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_R2.fastq.gz")
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "Bowtie", i)  +".sam"
    bowtie_align_cmd = ["bowtie", "-t", "-p", "18", "-x", reference_genome, "-1", input_fastq_r1, 
                         "-2", input_fastq_r2, "-S", output_sam]
    subprocess.run(bowtie_align_cmd, check=True)




#%% count features with htseq-count

try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "HTseq"))
except OSError as error:
    print(error)


bam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2")
bam_files = [bam for bam in os.listdir(bam_folder) if (bam.endswith(".bam") &  bam.startswith("sorted"))]  


for b in bam_files:  
    # Definir los archivos de entrada y salida
    bam_file =  os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2", b)
    output_htseq = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "HTseq", b)  +".counts"

    # Comando para HTSeq-count
    #command_htseq = ["htseq-count", "-f", "bam", "-t", 
     #          "gene", "-i", "gene_id", "-n", "20", "-r", "pos", bam_file, reference_genome_gtf]

    #SE data
    command_htseq = ["htseq-count", "-f", "bam", "-t", 
                     "gene", "-i", "gene_id", "-n", "20", bam_file, reference_genome_gtf]

    # Ejecutar el comando y redirigir la salida al archivo de conteo
    with open(output_htseq, "w") as output:
        subprocess.run(command_htseq, stdout=output, check=True)


#--num-threads  -n <n>, --nprocesses=<n>
#If the reads are PE you must specify 
#"-r", "name"
# "-s", "no" for non-stranded data

#%% use deseq2 to determine if the genes are differentially expressed

try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "DESeq2"))
except OSError as error:
    print(error)
DESeq2_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "DESeq2")


#Convert the data from HTseq to the format required by DESeq2
HTseq_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "HTseq")
HTseq_files = [HTseq for HTseq in os.listdir(HTseq_folder) if HTseq.endswith(".counts") ]  


dfs_HTseq = []
for HT in HTseq_files:  
    df = pd.read_table(os.path.join(HTseq_folder, HT), sep="\t", header=None)
    df = df.rename(columns={0: 'Gene_ID'}) 
    match1 = re.search(r"SC([^_]+)", HT)
    if match1:
        resultado = match1.group(0)
    df = df.rename(columns={1: resultado})
    df.set_index('Gene_ID', inplace=True)
    dfs_HTseq.append(df)

HT_concat = pd.concat(dfs_HTseq, axis = 1)

rows_eliminate = df.index[-5:]  # Get the index of the last 5 rows
HT_concat.drop(rows_eliminate, inplace=True)
HT_concat = HT_concat.transpose()
#Sort HTseq-counts results to coincide with the metadata observations

desired_index = ['SC24H-1', 'SC24H-2', 'SC48H-1', 'SC48H-2', 'SC72H-1', 'SC72H-2', 'SC96H-1', 'SC96H-2']
HT_concat= HT_concat.reindex(index=desired_index)

#save the results
HT_file = os.path.join(DESeq2_folder, "HTseq.txt")
HT_concat.to_csv(HT_file, sep = "\t", index=True)

#Open the metadata file
metadata_file = DESeq2_folder+"/metadata.txt"
metadata = pd.read_table(metadata_file, sep = "\t")
metadata.set_index('Sample', inplace=True)


#Filter out genes that have less than 5 read counts in total
genes_to_keep = HT_concat.columns[HT_concat.sum(axis=0) >= 10]
HT_concat = HT_concat[genes_to_keep]



# DeseqDataSet fits dispersion and log-fold change (LFC) parameters from
#the data, and stores them
dds = DeseqDataSet(
    counts=HT_concat,
    clinical=metadata,
    design_factors="Condition",
    refit_cooks=True,
    n_cpus=8,)

#Once a DeseqDataSet was initialized, we may run the deseq2() method to 
#fit dispersions and LFCs.
dds.deseq2()

#save the results
with open(os.path.join(DESeq2_folder, "dds_detailed_pipe.pkl"), "wb") as f:
        pkl.dump(dds, f)

#Now that dispersions and LFCs were fitted, we may proceed with
#statistical tests to compute p-values and adjusted p-values for 
#differential expresion. 
stat_res = DeseqStats(dds, n_cpus=8)

#PyDESeq2 computes p-values using Wald tests. 
#This can be done using the summary() method
stat_res.summary()

summary_DESeq2= stat_res.results_df

#For visualization or post-processing purposes, it might be suitable
#to perform LFC shrinkage. 
lfc_shrink = stat_res.lfc_shrink(coeff="Condition_Late-exponential_vs_Early-exponential")


#%% plots
## PCA
sc.tl.pca(dds)
sc.pl.pca(dds, color = 'Condition', size = 200)

#%% Volcano
FDR = summary_DESeq2['pvalue'] #cambiar por padjus
FDR = FDR.values
log_FDR = np.log10(FDR)

log2_FC = summary_DESeq2['log2FoldChange']
log2_FC = log2_FC.values
x_thr= 1
y_thr = -np.log10(0.3)

f, ax = plt.subplots()
ax.scatter(log2_FC[(log2_FC < -x_thr) & (-log_FDR > y_thr)], 
           -log_FDR[(log2_FC < -x_thr) & (-log_FDR > y_thr)], color='crimson', label='DE')

ax.scatter(log2_FC[(log2_FC > x_thr) & (-log_FDR > y_thr)], 
           -log_FDR[(log2_FC > x_thr) & (-log_FDR > y_thr)], color='dodgerblue', label='DE')

ax.scatter(log2_FC[(log2_FC < -x_thr) & (-log_FDR < y_thr)], 
           -log_FDR[(log2_FC < -x_thr) & (-log_FDR < y_thr)], color='darkgrey', label='not DE')

ax.scatter(log2_FC[(log2_FC > x_thr) & (-log_FDR < y_thr)], 
           -log_FDR[(log2_FC > x_thr) & (-log_FDR < y_thr)], color='darkgrey')

ax.scatter(log2_FC[(log2_FC > -x_thr) & (log2_FC < x_thr)], 
           -log_FDR[(log2_FC > -x_thr) & (log2_FC < x_thr)], 
           color='darkgrey')

ax.grid()
ax.set_xlabel("$log_2$ Fold Change")
ax.set_ylabel("$-log_{10}$ FDR")
ax.axvline(x=x_thr, color='green', linestyle='--')
ax.axvline(x=-x_thr, color='green', linestyle='--')
ax.axhline(y=y_thr, color='green', linestyle='--')
ax.legend()

# labels for the genes that have a FC > (1.5) and FC < -1.5
FC_thr = 4
summary_DESeq2_filtered = summary_DESeq2[(summary_DESeq2['log2FoldChange'] > FC_thr) | (summary_DESeq2['log2FoldChange'] < -FC_thr)]

labels = summary_DESeq2_filtered.index.tolist()
x_coords = list(summary_DESeq2_filtered['log2FoldChange'])
y_coords = list(summary_DESeq2_filtered['pvalue'])

for label, x, y in zip(labels, x_coords, y_coords):
    ax.annotate(label, (x, -np.log10(y)),textcoords="data", 
                xytext=(x, -np.log10(y)), ha='left', fontsize = 8)


#%% MA plot

log2_FC[(log2_FC < -x_thr) & (-log_FDR > y_thr)]

# Obtener los valores de basemean y log2 fold change
DE_genes = summary_DESeq2[summary_DESeq2['pvalue'] < 0.05]
non_DE = summary_DESeq2[summary_DESeq2['pvalue'] >= 0.05]

basemean_DE = DE_genes.baseMean
basemean_non_DE = non_DE.baseMean

log2_FC_DE = DE_genes.log2FoldChange
log2_FC_non_DE = non_DE.log2FoldChange


# Crear la figura y el eje
fig2, ax = plt.subplots()
ax.scatter(basemean_DE, log2_FC_DE, color='red', label = "DE")
ax.scatter(basemean_non_DE, log2_FC_non_DE, color='black', label = "non DE")
ax.legend()
ax.axhline(y=0, color='red', linestyle='--')
ax.grid()
ax.set_xlabel('mean of normalized counts')
ax.set_ylabel('$log_2$ Fold Change')



#%% Align sequences to the chromosome and plasmid individually with with  Bowtie2

# generate the index genome file
# Comando Bowtie2 para crear la base de datos de referencia

#CHR
bowtie2_build_cmd = ["bowtie2-build", reference_chr_fasta, reference_chr]
subprocess.run(bowtie2_build_cmd, check=True)

#PLASMID
bowtie2_build_cmd = ["bowtie2-build", reference_plasmid_fasta, reference_plasmid]
subprocess.run(bowtie2_build_cmd, check=True)


fastq_filtered_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp")
fastq_filtered = [fastq for fastq in os.listdir(fastq_filtered_folder) if (fastq.endswith(".gz") or fastq.endswith(".fastq")) ]  






# Alignment for PE data CHR
try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_CHR"))
except OSError as error:
    print(error)

for i in experiments:  
    # Comamand for 
    input_fastq_r1 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_R1.fastq.gz")
    input_fastq_r2 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_R2.fastq.gz")
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "Bowtie2_CHR", i)  +".sam"
    bowtie2_align_cmd = ["bowtie2", "-p", "18", "-x", reference_chr, "-1", input_fastq_r1, 
                         "-2", input_fastq_r2, "-S", output_sam]
    subprocess.run(bowtie2_align_cmd, check=True)
    
    
# Alignment for PE data PLASMID
try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_PLASMID"))
except OSError as error:
    print(error)

for i in experiments:  
    # Comamand for 
    input_fastq_r1 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_R1.fastq.gz")
    input_fastq_r2 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_R2.fastq.gz")
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "Bowtie2_PLASMID", i)  +".sam"
    bowtie2_align_cmd = ["bowtie2", "-p", "18", "-x", reference_plasmid, "-1", input_fastq_r1, 
                         "-2", input_fastq_r2, "-S", output_sam]
    subprocess.run(bowtie2_align_cmd, check=True)
    
#%% Convert SAM files to BAM files 
  
# Alignment for PE CHR
alignment_sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_CHR")
alignment_sam = [sam for sam in os.listdir(alignment_sam_folder) if sam.endswith(".sam") ]  
for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    output_bam = os.path.join(alignment_sam_folder, os.path.splitext(alignment)[0]+".bam")
    convert_sam_to_bam(input_sam, output_bam)

  
# Alignment for PE plasmid
alignment_sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_PLASMID")
alignment_sam = [sam for sam in os.listdir(alignment_sam_folder) if sam.endswith(".sam") ]  
for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    output_bam = os.path.join(alignment_sam_folder, os.path.splitext(alignment)[0]+".bam")
    convert_sam_to_bam(input_sam, output_bam)
#%% Convert unsorted to sorted BAM files

sorted_bam(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_CHR"))
sorted_bam(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_PLASMID"))

#%%  Delete sam files

delete_files_aln(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_CHR"), ".sam")
delete_files_aln(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_PLASMID"), ".sam")

#%% generate .bai files


#CHR
alignment_sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_CHR")
alignment_bam_sorted = [bam for bam in os.listdir(alignment_sam_folder) if bam.startswith("sorted")]  
bam_index_file(alignment_sam_folder, alignment_bam_sorted)


#PLASMID
alignment_sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_PLASMID")
alignment_bam_sorted = [bam for bam in os.listdir(alignment_sam_folder) if bam.startswith("sorted")]  
bam_index_file(alignment_sam_folder, alignment_bam_sorted)


#%%

from IPython import get_ipython
get_ipython().magic('reset -sf')


import subprocess
import os
import shutil
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt


apero_folder=  os.path.join(cur_dir, "results_rna_seq", my_experiment, "APERO")
plasmid_files = [file for file in os.listdir(apero_folder) if (file.endswith(".csv") &  file.startswith("PLASMID"))]  
chr_files = [file for file in os.listdir(apero_folder) if (file.endswith(".csv") &  file.startswith("CHR"))]  


#Open results from replicate 1
for i in plasmid_files
    df1 = pd.read_table(os.path.join(apero_folder, i), sep=",")
    #Filter for candidates that start with more than 10 reads
    df1_coverage = df1[df1['freq'] >= 10]
    #Discard internall small transcripts located inside CDSs
    df1_internal = df1_coverage[df1_coverage['Class'] != 'Ai']





# Assume you have a list of file names in a folder
folder = 'your_folder_with_files'
files = [file for file in os.listdir(folder) if file.endswith('.csv')]
full_paths = [os.path.join(folder, file) for file in files]

# Initialize the DataFrame result_final
result_final = pd.DataFrame()

# Iterate over each file and create a DataFrame for each one
for file_path in full_paths:
    df = pd.read_csv(file_path)

    # Perform the comparison and information preservation procedure
    for index, reference_row in df.iterrows():
        reference_position = reference_row['Position']
        lg_maximum = reference_row['lg']

        matching_rows = []
        for other_df in dataframes:
            matching_rows.extend(other_df[
                (other_df['Position'] >= (reference_position - 10)) &
                (other_df['Position'] <= (reference_position + 10))
            ])

        if matching_rows:
            df_maximum = max(matching_rows, key=lambda x: x['lg'])
            result_final = pd.concat([result_final, df_maximum])

# Reset the indices of the final DataFrame
result_final.reset_index(drop=True, inplace=True)




    #Open results from replicate 1
    df2 = pd.read_table(os.path.join(apero_folder, "PLASMID_sorted_trimmedSC24H-2_S399_L001.csv"),
                   sep=",")
    #Filter for candidates that start with more than 10 reads
    df2_coverage = df2[df2['freq'] >= 10]
    #Discard internall small transcripts located inside CDSs
    df2_internal = df2_coverage[df2_coverage['Class'] != 'Ai']
    

from BCBio import GFF
from Bio import SeqIO

#Extract fasta sequence of each replicon
reference_genome_gbk = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/GCF_005519465.1_ASM551946v1_genomic.gbff"



for seq_req in SeqIO.parse(reference_genome_gbk , "genbank"):
    print(seq_req.id)
    SeqIO.write(seq_req, "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/"+ seq_req.id + ".gbk", "genbank")
    in_file = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/"+ seq_req.id+ ".gbk"
    out_file = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/"+ seq_req.id+ ".gff3"
    in_handle = open(in_file)
    out_handle = open(out_file, "w")

    GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)

    in_handle.close()
    out_handle.close()


