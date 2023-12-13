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
reference_genome = "/home/usuario/Documentos/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/GCF_005519465.1_ASM551946v1_genomic.fna"
cur_dir = os.getcwd()
#path_rawdata = os.path.join(cur_dir, "raw_data")
path_rawdata = os.path.join(cur_dir, "raw_data", "RNA-seq_OMEGA")


fastq_files = [file for file in os.listdir(path_rawdata) if file.endswith(".gz") ]
fastq_files_path = []
for i in fastq_files:
    fastq_files_path.append(os.path.join(path_rawdata, i))



#%% Initial quality control of the data with fastqc

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
                   "--cut_right", "-j", j+".json", "-h", j+".html", "-w", "16"]

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

#%% Fastqc after rRNA removal

fastq_files_rRNA_removed= [file for file in os.listdir(fastq_rRNA_removed) if file.endswith(".gz") ]


fastQC_command_rRNA_removed =[]  
for i, j  in enumerate(fastq_files_rRNA_removed):
    fastQC_command_rRNA_removed.append(os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed", j))
fastQC_rRNA_removed=os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastQC_rRNA_removed")
try:
    os.mkdir(fastQC_rRNA_removed)
except OSError as error:
    print(error) 

   
fastQC_command_rRNA_removed.insert(0, "fastqc")
fastQC_command_rRNA_removed.append("-o")
fastQC_command_rRNA_removed.append(fastQC_rRNA_removed)
fastQC_command_rRNA_removed.append("-t")
fastQC_command_rRNA_removed.append("20")   #change the number of threads accordingly
#--noextract
subprocess.call(fastQC_command_rRNA_removed)


#Summarize all the fastqc results with multiQC
MultiQC(fastQC_rRNA_removed, "multiqc_fastp_rRNA_removal")


#%% reformat the output of sortmeRNA to have again the forward and reverse fastq files
#as indepedent files


for i, j  in enumerate(fastq_files_rRNA_removed):
   reformat_command = ["reformat.sh", 
                        "in="+os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed", j), 
                        "out1="+os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed", j.split(".fq")[0]+"_1.fastq.gz"),
                        "out2="+os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed", j.split(".fq")[0]+"_2.fastq.gz")]
   subprocess.call(reformat_command)



#%% Removeinterleaved files

character = '_'
for filename in os.listdir(fastq_rRNA_removed):
    if character not in filename:
        file_path = os.path.join(fastq_rRNA_removed, filename)
        # Remover el archivo
        os.remove(file_path)


#%% Align sequences to the reference genome  Bowtie2

try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2"))
except OSError as error:
    print(error)

fastq_filtered_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed")
fastq_filtered = [fastq for fastq in os.listdir(fastq_filtered_folder) if fastq.endswith(".gz") ]  

# generate the index genome file
# Comando Bowtie2 para crear la base de datos de referencia
bowtie2_build_cmd = ["bowtie2-build", reference_genome, "genome"]
subprocess.run(bowtie2_build_cmd, check=True)

#Drop the files that belongs to pseudomonas
experiments_Streptomyces = experiments
experiments_Streptomyces.remove('PACP72H-1')
experiments_Streptomyces.remove('PACP72H-2')

#fastq_filtered.remove('PACP72H-2_1.fastq.gz')
#fastq_filtered.remove('PACP72H-2_2.fastq.gz')

#SE data
for i in fastq_filtered:  
    fastq_file = os.path.join(cur_dir, "results_rna_seq", my_experiment, "cutadapt", i)
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "Bowtie2_cutadapt", i)  +".sam"
    bowtie2_align_cmd = ["bowtie2", "-p", "20", "-x",  "genome", "-U", fastq_file,  "-S",  output_sam]
    subprocess.run(bowtie2_align_cmd, check=True)
    
#PE data
for i in experiments_Streptomyces:  
    # Comamand for 
    input_fastq_r1 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed", i+"_1.fastq.gz")
    input_fastq_r2 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed", i+"_2.fastq.gz")
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "Bowtie2", i)  +".sam"
    bowtie2_align_cmd = ["bowtie2", "-p", "18", "-x", "genome", "-1", input_fastq_r1, 
                         "-2", input_fastq_r2, "-S", output_sam]
    subprocess.run(bowtie2_align_cmd, check=True)


#%% Align sequences to the reference genome  Bowtie2 without rRNA-removal

try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_fastp"))
except OSError as error:
    print(error)

fastq_filtered_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp")
fastq_filtered = [fastq for fastq in os.listdir(fastq_filtered_folder) if fastq.endswith(".gz") ]  

# generate the index genome file
# Comando Bowtie2 para crear la base de datos de referencia
bowtie2_build_cmd = ["bowtie2-build", reference_genome, "genome"]
subprocess.run(bowtie2_build_cmd, check=True)


#SE data
for i in fastq_filtered:  
    fastq_file = os.path.join(cur_dir, "results_rna_seq", my_experiment, "cutadapt", i)
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "Bowtie2_cutadapt", i)  +".sam"
    bowtie2_align_cmd = ["bowtie2", "-p", "20", "-x",  "genome", "-U", fastq_file,  "-S",  output_sam]
    subprocess.run(bowtie2_align_cmd, check=True)
    
#PE data
for i in experiments:  
    # Comamand for 
    input_fastq_r1 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_R1.fastq.gz")
    input_fastq_r2 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_R2.fastq.gz")
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "Bowtie2_fastp", i)  +".sam"
    bowtie2_align_cmd = ["bowtie2", "-p", "18", "-x", "genome", "-1", input_fastq_r1, 
                         "-2", input_fastq_r2, "-S", output_sam]
    subprocess.run(bowtie2_align_cmd, check=True)

#%% Convert SAM files to BAM files 

alignment_sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_fastp")
alignment_sam = [sam for sam in os.listdir(alignment_sam_folder) if sam.endswith(".sam") ]  


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

for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    output_bam = os.path.join(alignment_sam_folder, os.path.splitext(alignment)[0]+".bam")
    convert_sam_to_bam(input_sam, output_bam)

#%% Convert unsorted to sorted BAM files
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
