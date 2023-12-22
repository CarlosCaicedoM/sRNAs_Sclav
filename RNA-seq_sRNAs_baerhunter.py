#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 10:15:11 2023

@author: usuario
"""

#conda install -c bioconda fastqc
#sudo apt install sra-toolkit
#conda install -c bioconda -c conda-forge multiqc
#conda create -n cutadaptenv cutadapt
#sudo apt install cutadapt
# conda install -c bioconda fastp
#sudo apt install bwa
#sudo apt install samtools
#pip install HTSeq
#pip install pydeseq2
#conda install -c conda-forge scanpy python-igraph leidenalg
#pip install gseapy
#pip install goatools
#pip install pysam
#sudo apt install bowtie2

#%%

from IPython import get_ipython
get_ipython().magic('reset -sf')


import subprocess
import os
import shutil
import multiqc
import pandas as pd
import re
import pickle as pkl
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import scanpy as sc
import gseapy as gp
from gseapy.plot import gseaplot
import numpy as np
import matplotlib.pyplot as plt
import pysam
import gzip


#Create file to download your data

my_experiment = "MTB_baerhunter"
my_file = "/SraAccList_MTB_baerhunter.txt"
reference_genome_fna = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/MTB_annotations/GCF_000195955.2_ASM19595v2_genomic.fna"
reference_genome_gtf = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/MTB_annotations/GCF_000195955.2_ASM19595v2_genomic.gtf"
reference_genome = "MTB-H37Rv"



cur_dir = os.getcwd()
path_rawdata = os.path.join(cur_dir, "raw_data", my_experiment) 


try: 
    os.mkdir(path_rawdata) 
except OSError as error: 
    print(error) 

my_data = cur_dir+my_file
#If you want to remove all whitespace characters (newlines and spaces) from the end of each line
with open(my_data) as f:
    lines = [line.rstrip() for line in f]
Accessions  = lines[0:-1]

# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in Accessions:
    print ("Currently downloading: " + sra_id)
    fasterq = "fasterq-dump " + sra_id + " -O " + path_rawdata
    print ("The command used was: " + fasterq)
    subprocess.call(fasterq, shell=True)
	
fastq_files = [file for file in os.listdir(path_rawdata) if file.endswith(".fastq") ]
fastqc_outdir=os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastqc")

try:
    os.mkdir("results_rna_seq")
except OSError as error:
    print(error)  

try:
    os.mkdir(os.path.join("results_rna_seq", my_experiment))
except OSError as error:
    print(error)

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
                   "-i", path_rawdata +"/"+j, "-o",  
                   os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+j), 
                   "-j", j+".json", "-h", j+".html", "-w", "16", "--cut_right"]

    subprocess.call(fastp_command)
"""
#If paired-end

for i,j  in enumerate(Accessions):
    fastp_command=["fastp", 
                   "-i", path_rawdata+"/"+j+"_1.fastq",
                   "-I",  path_rawdata+"/"+j+"_2.fastq",               
                   "-o", os.path.join(cur_dir, "results_rna_seq",
                                my_experiment, "fastp", "trimmed_"+j+"_1.fastq"), 
                   "-O", os.path.join(cur_dir, "results_rna_seq",
                                my_experiment, "fastp", "trimmed_"+j+"_2.fastq"),
                   "-j", j+".json", "-h", j+".html", "-w", "16", "--cut_right", 
                   "--detect_adapter_for_pe"]



    subprocess.call(fastp_command)
"""    
    
#Move the reports to the fastp folder
reports = [rep for rep in os.listdir(os.curdir) if (rep.endswith(".html") or rep.endswith(".json"))]  
for file_name in reports:
    shutil.move(file_name, os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp"))


#Move the reports to the fastp folder
reports = [rep for rep in os.listdir(os.curdir) if (rep.endswith(".html") or rep.endswith(".json"))]  
for file_name in reports:
    shutil.move(file_name, os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp"))


"""fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
for paired end data modify accordingly"""

#adapter trimming is enabled by default
#quality filtering is enabled by default. 
 #length filtering is enabled by default.
#-w number of threads 
#-detect_adapter_for_pe        
#by default, the auto-detection for adapter is for SE data input only,
#turn on this option to enable it for PE data
#WARNING: fastp uses up to 16 threads 

"""
-r, --cut_right move a sliding window from front to tail, 
if meet one window with mean quality < threshold, drop the bases in the 
window and the right part, and then stop. Use cut_right_window_size to set 
the widnow size, and cut_right_mean_quality to set the mean quality threshold. 
This is similar as the Trimmomatic SLIDINGWINDOW method

 "--trim_poly_g"
"""
#%% Fastqc after trimming

fastq_files_trimmed = [file for file in os.listdir(fastp_dir) if (file.endswith(".gz") or file.endswith(".fastq"))]
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

#%% Align sequences to the reference genome  Bowtie2

# generate the index genome file
# Comando Bowtie2 para crear la base de datos de referencia

bowtie2_build_cmd = ["bowtie2-build", reference_genome_fna, reference_genome]
subprocess.run(bowtie2_build_cmd, check=True)


# Alignment for SE data

try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2"))
except OSError as error:
    print(error)

fastq_filtered_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp")
fastq_filtered = fastq_files_trimmed 

for i in fastq_filtered:  
    fastq_file = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", i)
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "Bowtie2", i)  +".sam"
    bowtie2_align_cmd = ["bowtie2", "-p", "18", "-x",  reference_genome, "-U", fastq_file,  "-S",  output_sam]
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
    

alignment_sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2")
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


for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    is_sorted = check_bam_order(input_sam)

if is_sorted:
    print("The BAM file is already sorted.")
else:
    sorted_bam(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2"))


#%%  Delete sam files and unsorted bam


def delete_files_aln(folder, suffix):
    files = os.listdir(folder)
    for file in files:
        path_file = os.path.join(folder, file)
        if file.endswith(suffix) or (file.endswith(".bam") and not file.startswith("sorted")):
            os.remove(path_file)

delete_files_aln(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2"), ".sam")

#%% Align sequences to the reference genome with BWA

try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "BWA"))
except OSError as error:
    print(error)


# generate the index genome file
command = ['bwa', 'index', reference_genome_fna]
subprocess.run(command, check=True)


#SE data
for i in fastq_filtered:  
    # Comamand for BWA
    fastq_file = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", i)
    command = ["bwa", "mem", reference_genome_fna, fastq_file]
    # redirect the output to the SAM file
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "BWA", i)  +".sam"
    
    with open(output_sam, "w") as output_file:
        subprocess.run(command, stdout=output_file, check=True)
    


#PE data

for i in Accessions:  
    # Comamand for BWA
    fastq_file_1 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_1.fastq")
    fastq_file_2 = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+i+"_2.fastq")
    command = ["bwa", "mem", reference_genome, fastq_file_1, fastq_file_2]
    # redirect the output to the SAM file
    output_sam = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "BWA", i)  +".sam"
    
    with open(output_sam, "w") as output_file:
        subprocess.run(command, stdout=output_file, check=True)
        
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
