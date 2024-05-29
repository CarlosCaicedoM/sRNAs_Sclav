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

from pydeseq2.default_inference import DefaultInference


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
reference_genome_bed = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/GCF_005519465.1_ASM551946v1_genomic.bed"
reference_genome = "SCLAV"

reference_chr_fasta = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/NZ_CP027858.1/NZ_CP027858.1.fasta"
reference_plasmid_fasta = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/NZ_CP027859.1/NZ_CP027859.1.fasta"
reference_chr = "SCLAV_CHR"
reference_plasmid = "SCLAV_PLASMID"

#%% main
cur_dir = os.getcwd()
#path_rawdata = os.path.join(cur_dir, "raw_data")
path_rawdata = os.path.join(cur_dir, "raw_data", "RNA-seq_OMEGA")


fastq_files = [file for file in os.listdir(path_rawdata) if file.endswith(".gz") ]
fastq_files_path = []
for i in fastq_files:
    fastq_files_path.append(os.path.join(path_rawdata, i))

#Initial quality control of the data with fastqc

try:
    os.mkdir("results_rna_seq")
except OSError as error:
    print(error)  

try:
    os.mkdir(os.path.join("results_rna_seq", my_experiment))
except OSError as error:
    print(error)

#Initial quality control with FASTQC
for i in fastq_files:
    FASTQC_conda_env("FASTQC", path_rawdata, "fastqc", 18)

#Summarize the initial qulity control results with MultiQC
MultiQC(fastqc_outdir, "multiqc")

# Quality control and filtering by fastp
fastp_PE(fastq_files, fastq_files_path, "fastp", 16, "_")

# Fastqc after trimming
fastp_dir =os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "fastp")
fastq_files_trimmed = [file for file in os.listdir(fastp_dir) if (file.endswith(".gz") or file.endswith(".fastq")) ]
fastq_files_trimmed_path = []
for i in fastq_files_trimmed_path:
    fastq_files_path.append(os.path.join(fastp_dir, i))
for i in fastq_files_trimmed:
    FASTQC_conda_env("FASTQC", fastq_files_trimmed_path, "fastQC_trimmed", 18)

#Summarize all the fastqc results with multiQC
fastQC_trimmed = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "fastQC_trimmed")
MultiQC(fastQC_trimmed, "multiqc_fastp")


#%  Remove ribosomal RNA sequences
sortmeRNA_PE(fastq_files_trimmed, "fastp", "fastq_rRNA_removed", 18, "_")

#reformat the output of sortmeRNA to have again the forward and reverse fastq files
#as indepedent files
fastq_rRNA_removed=os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "fastq_rRNA_removed")
fastq_files_rRNA_removed= [file for file in os.listdir(fastq_rRNA_removed) if file.endswith(".gz") ]
reformat_BBmap(fastq_rRNA_removed, "fastq_rRNA_removed")

# Removeinterleaved files
remove_interleaved_files(fastq_rRNA_removed)

#Fastqc after rRNA removal
fastq_files_RNA_removed_path = []
for i in fastq_files_rRNA_removed:
    fastq_files_path.append(os.path.join(fastq_rRNA_removed, i))
for i in fastq_files_RNA_removed_path:
    FASTQC_conda_env("FASTQC", fastq_files_RNA_removed_path, "fastQC_rRNA_removed", 18)

#Summarize all the fastqc results with multiQC
fastQC_rRNA_removed = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "fastQC_rRNA_removed")
MultiQC(fastQC_rRNA_removed, "multiqc_fastp_rRNA_removal")


# Align sequences to the reference genome  Bowtie2
## generate the index genome file
generate_index_genome_bowtie2(reference_genome_fna, reference_genome)
 
## PE_DATA without RNA removal
fastq_filtered_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp")
pattern = re.compile(r".*_[12]\.(?:fastq\.gz|gz|fastq)$")
fastq_filtered = [file for file in os.listdir(fastq_filtered_folder) if pattern.match(file)]
align_bowtie2_PE(fastq_filtered, fastq_filtered_folder, "Bowtie2", reference_chr, "_", 18)

## Alignment for PE data WITH rRNA removal 
fastq_filtered_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastq_rRNA_removed")
pattern = re.compile(r".*_[12]\.(?:fastq\.gz|gz|fastq)$")
fastq_filtered = [file for file in os.listdir(fastq_filtered_folder) if pattern.match(file)]
align_bowtie2_PE(fastq_filtered, fastq_filtered_folder, "Bowtie2_rRNA_removed", reference_chr, "_", 18)

# Convert SAM files to BAM files 
    
## Alignment for PE data WITHOUT rRNA removal 
alignment_sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2")
alignment_sam = [sam for sam in os.listdir(alignment_sam_folder) if sam.endswith(".sam") ]  
for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    output_bam = os.path.join(alignment_sam_folder, os.path.splitext(alignment)[0]+".bam")
    convert_sam_to_bam(input_sam, output_bam)

## Alignment for PE data WITH rRNA removal 
alignment_sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_rRNA_removed")
alignment_sam = [sam for sam in os.listdir(alignment_sam_folder) if sam.endswith(".sam") ]  
for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    output_bam = os.path.join(alignment_sam_folder, os.path.splitext(alignment)[0]+".bam")
    convert_sam_to_bam(input_sam, output_bam)
    
#Convert unsorted to sorted BAM files
           
sorted_bam(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2"))
sorted_bam(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_rRNA_removed"))

#Check if the rRNA was succesfully removed by counting the number of 
#reads aligned to rRNA genes
alignment_bam_sorted = [bam for bam in os.listdir(alignment_sam_folder) if bam.startswith("sorted")]  
convert_gff_to_bed(reference_genome_gtf, reference_genome_bed)
reference_genome_bed_rRNA = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/GCF_005519465.1_ASM551946v1_rRNA.bed"

calculate_rrna_percentage(alignment_bam_sorted, alignment_sam_folder, reference_genome_bed_rRNA)

#access a BAM index file
bam_index_file(alignment_sam_folder, alignment_bam_sorted)


#Delete sam files

delete_files_aln(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2"), ".sam")
delete_files_aln(os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2_rRNA_removed"), ".sam")


# Align sequences to the reference genome  BWA
## generate the index genome file
generate_index_genome_BWA(reference_genome_fna)

fastq_filtered_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp")
pattern = re.compile(r".*_[12]\.(?:fastq\.gz|gz|fastq)$")
fastq_filtered = [file for file in os.listdir(fastq_filtered_folder) if pattern.match(file)]
align_BWA_PE(fastq_filtered, fastq_filtered_folder, "BWA", reference_genome_fna , "_", 18)

#Convert SAM files to BAM files 
alignment_sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "BWA")
alignment_sam = [sam for sam in os.listdir(alignment_sam_folder) if sam.endswith(".sam") ]  
for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    output_bam = os.path.join(alignment_sam_folder, os.path.splitext(alignment)[0]+".bam")
    convert_sam_to_bam(input_sam, output_bam)


#Convert unsorted to sorted BAM files
sorted_bam(os.path.join(cur_dir, "results_rna_seq", my_experiment, "BWA"))


#Delete sam files and unsorted bam
delete_files_aln(os.path.join(cur_dir, "results_rna_seq", my_experiment, "BWA"), ".sam")

#access a BAM index file
alignment_sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "BWA")
alignment_bam_sorted = [bam for bam in os.listdir(alignment_sam_folder) if bam.startswith("sorted")]  
bam_index_file(alignment_sam_folder, alignment_bam_sorted)

# Align sequences to the reference genome  Bowtie 

# generate the index genome file
generate_index_genome_bowtie(reference_genome_fna, reference_genome)
    
align_bowtie_PE(fastq_filtered, fastq_filtered_folder, "Bowtie", reference_genome_fna, "_", 18)

    
# count features with htseq-count


 


#%% Functions
        
def FASTQC_conda_env(conda_env_name, fastq_files_path, results_folder, threads):
    fastqc_outdir=os.path.join(os.getcwd(), "results_rna_seq", my_experiment, results_folder)
    try:
        os.mkdir(fastqc_outdir)
    except OSError as error:
        print(error)        
    fastQC_command =[]  
    for i in fastq_files:
        fastQC_command.append(i)
    fastQC_command.insert(0, "fastqc")
    fastQC_command.append("-o")
    fastQC_command.append(fastqc_outdir)
    fastQC_command.append("-t")
    fastQC_command.append(str(threads))
    # Use conda run to execute the command in the specified environment
    full_command = ["conda", "run", "-n", conda_env_name] + fastQC_command
    subprocess.run(full_command)

def FASTQC(fastq_files_path, results_folder, threads):
    fastqc_outdir=os.path.join(os.getcwd(), "results_rna_seq", my_experiment, results_folder)
    try:
        os.mkdir(fastqc_outdir)
    except OSError as error:
        print(error)
    fastQC_command =[]  
    for i in fastq_files:
        fastQC_command.append(i)
    fastQC_command.insert(0, "fastqc")
    fastQC_command.append("-o")
    fastQC_command.append(fastqc_outdir)
    fastQC_command.append("-t")
    fastQC_command.append(str(threads))   #change the number of threads accordingly
    #--noextract
    subprocess.call(fastQC_command)

def MultiQC(inputdir, outdir):
    try:
        os.mkdir(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, outdir))
    except OSError as error:
        print(error)
        
    multiQC_command = ["multiqc", inputdir]
    multiQC_command.append("-n")
    multiQC_command.append(my_experiment)
    multiQC_command.append("-o")
    multiQC_command.append(os.path.join(os.getcwd(), "results_rna_seq", my_experiment,  outdir))
    multiQC_command.append("-f")
    multiQC_command.append("--profile-runtime") 
    subprocess.call(multiQC_command)


def fastp_SE(fastq_files, fastq_files_path, results_folder, threads):
    fastp_dir =os.path.join(os.getcwd(), "results_rna_seq", my_experiment, results_folder)
    try:
      os.mkdir(fastp_dir)
    except OSError as error:
        print(error)    
    for i,j  in enumerate(fastq_files):
        fastp_command=["fastp", 
                       "-i", fastq_files_path[i], "-o",  
                       os.path.join(os.getcwd(), "results_rna_seq", my_experiment, results_folder, "trimmed_"+j),
                       "-j", j+".json", "-h", j+".html", "-w", str(threads), 
                       "--cut_right", "--trim_poly_g", "--trim_poly_x"]
    
        subprocess.call(fastp_command)
    #Move the reports to the fastp folder
    reports = [rep for rep in os.listdir(os.curdir) if (rep.endswith(".html") or rep.endswith(".json"))]  
    for file_name in reports:
        shutil.move(file_name, os.path.join(cur_dir, "results_rna_seq", my_experiment, results_folder))
        

def fastp_PE(fastq_files, fastq_files_path, results_folder, threads, end):
    fastp_dir =os.path.join(os.getcwd(), "results_rna_seq", my_experiment, results_folder)
    try:
      os.mkdir(fastp_dir)
    except OSError as error:
        print(error)   
        
    experiments = []
    for file in fastq_files:
        parts = file.split(end)
        name = parts[0]
        experiments.append(name)
    experiments = set(experiments) 
    
    for i,j  in enumerate(experiments):
        input_1 = next(file for file in fastq_files_path if re.search(f"{j}{end}1\.(fastq|gz|fastq\.gz)$", file))
        input_2 = next(file for file in fastq_files_path if re.search(f"{j}{end}2\.(fastq|gz|fastq\.gz)$", file))
        output_1 = os.path.join(os.getcwd(), "results_rna_seq",
                     my_experiment, results_folder, "trimmed-"+j+"_1.fastq.gz")
        output_2 = os.path.join(os.getcwd(), "results_rna_seq",
                     my_experiment, results_folder, "trimmed-"+j+"_2.fastq.gz")
        fastp_command=["fastp", 
                       "-i", input_1,
                       "-I", input_2,                                           
                       "-o", output_1, 
                       "-O", output_2,
                       "-j", j+".json", "-h", j+".html", "-w", str(threads), 
                       "--cut_right", "--trim_poly_g", "--trim_poly_x"]
    
        subprocess.call(fastp_command)
    #Move the reports to the fastp folder
    reports = [rep for rep in os.listdir(os.curdir) if (rep.endswith(".html") or rep.endswith(".json"))]  
    for file_name in reports:
        shutil.move(file_name, os.path.join(os.getcwd(), "results_rna_seq", my_experiment, results_folder))

def sortmeRNA_PE(fastq_files, fastq_files_folder, results_folder, threads, end):
    fastq_rRNA_removed=os.path.join(os.getcwd(), "results_rna_seq", my_experiment, results_folder)
    try:
        os.mkdir(fastq_rRNA_removed)
    except OSError as error:
        print(error)   
        
    experiments = []
    for file in fastq_files:
        parts = file.split(end)
        name = parts[0]
        experiments.append(name)
    experiments = set(experiments) 
    
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
                             "--reads", os.path.join(os.getcwd(), "results_rna_seq",
                                          my_experiment, fastq_files_folder, "trimmed-"+j+"_1.fastq.gz"), 
                             "--reads", os.path.join(os.getcwd(), "results_rna_seq",
                                          my_experiment, fastq_files_folder, "trimmed-"+j+"_2.fastq.gz"),
                             "--other",  os.path.join(fastq_rRNA_removed, j), 
                             "--fastx", "--threads",  str(threads), "--paired_in", 
                             "--workdir", j] 
    
        subprocess.call(sortmeRNA_command)
        # Move the directory to fastq_rRNA_removed
        shutil.move(j, fastq_rRNA_removed)
        
def reformat_BBmap(fastq_files_rRNA_removed, folder):
    for i, j  in enumerate(fastq_files_rRNA_removed):
       reformat_command = ["reformat.sh", 
                            "in="+os.path.join(os.getcwd(), "results_rna_seq", my_experiment, folder, j), 
                            "out1="+os.path.join(os.getcwd(), "results_rna_seq", my_experiment, folder, j.split(".fq")[0]+"_1.fastq.gz"),
                            "out2="+os.path.join(os.getcwd(), "results_rna_seq", my_experiment, folder, j.split(".fq")[0]+"_2.fastq.gz")]
       subprocess.call(reformat_command)

def remove_interleaved_files(fastq_rRNA_removed):
    character = '.fastq.'
    for filename in os.listdir(fastq_rRNA_removed):
        if character not in filename:
            file_path = os.path.join(fastq_rRNA_removed, filename)
            # Remover el archivo
            os.remove(file_path)

def generate_index_genome_bowtie2(reference_fasta, reference_genome):
    bowtie2_build_cmd = ["bowtie2-build", reference_fasta, reference_genome]
    subprocess.run(bowtie2_build_cmd, check=True)

def align_bowtie2_SE(fastq_files, path_rawdata, folder_results, reference_genome, threads):
    try:
      os.mkdir(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, folder_results))
    except OSError as error:
        print(error)
    for i in fastq_files:  
        fastq_file = os.path.join(path_rawdata, i)
        output_sam = os.path.join(os.getcwd(), "results_rna_seq",
                                  my_experiment, folder_results, i)  +".sam"
        bowtie2_align_cmd = ["bowtie2", "-p", str(threads),
                             "-x",  reference_genome,
                             "-U", fastq_file,
                             "-S",  output_sam]
        subprocess.run(bowtie2_align_cmd, check=True)


def align_bowtie2_PE(fastq_files, path_rawdata, folder_results, reference_genome, end, threads):
    try:
      os.mkdir(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, folder_results))
    except OSError as error:
        print(error)

    experiments = []
    for file in fastq_files:
        parts = file.split(end)
        name = parts[0]
        experiments.append(name)
    
    experiments = set(experiments) 
    
    for i in experiments:  
       # Comamand for 
       input_fastq_r1 = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, path_rawdata, i+"_1.fastq.gz")
       input_fastq_r2 = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, path_rawdata, i+"_2.fastq.gz")
       output_sam = os.path.join(os.getcwd(), "results_rna_seq",
                                 my_experiment, folder_results, i)  +".sam"
       bowtie2_align_cmd = ["bowtie2", "-p", str(threads), "-x", reference_genome, "-1", input_fastq_r1, 
                            "-2", input_fastq_r2, "-S", output_sam]
       subprocess.run(bowtie2_align_cmd, check=True)


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


def convert_gff_to_bed(gff_file, bed_file):
    cmd = f"gtf2bed < {gff_file} > {bed_file}"
    subprocess.run(cmd, shell=True)

def bed_rna(reference_genome_bed, reference_genome_bed_rRNA):
    command = f"grep 'gbkey \"rRNA\"' {reference_genome_bed} > {reference_genome_bed_rRNA}"
    subprocess.run(command, shell=True)

def bam_index_file(alignment_sam_folder, alignment_bam_sorted):
    for i in alignment_bam_sorted:
        bam_file = os.path.join(alignment_sam_folder, i)
        command_index = ["samtools", "index", bam_file]    
        subprocess.run(command_index)

def calculate_rrna_percentage(alignment_bam_sorted, alignment_sam_folder, reference_genome_bed_rRNA):
    for i in alignment_bam_sorted:
        bam_file = pysam.AlignmentFile(os.path.join(alignment_sam_folder, i), 'rb')
        bed_file = pybedtools.BedTool(reference_genome_bed_rRNA)
        
        # Filter reads aligning to rRNA regions
        filtered_reads = bed_file.intersect(os.path.join(alignment_sam_folder, i), u=True)
        
        # Count reads in BAM files
        total_reads = sum(1 for read in bam_file)
        rRNA_reads = sum(1 for read in filtered_reads)
        
        # Calculate percentage
        percentage = 100.0 * rRNA_reads / total_reads
        
        print(f'The percentage of reads aligned to rRNA is: {percentage}%')
        
        # Close BAM file
        bam_file.close()


def delete_files_aln(folder, suffix):
    files = os.listdir(folder)
    for file in files:
        path_file = os.path.join(folder, file)
        if file.endswith(suffix) or (file.endswith(".bam") and not file.startswith("sorted")):
            os.remove(path_file)

def generate_index_genome_BWA(reference_genome_fasta):
    command = ['bwa', 'index', reference_genome_fasta]
    subprocess.run(command, check=True)


def align_BWA_SE(fastq_files, path_rawdata, folder_results, reference_genome_fasta, threads):
    try:
      os.mkdir(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, folder_results))
    except OSError as error:
        print(error)

    for i in fastq_files:  
        # Comamand for BWA
        fastq_file = os.path.join(path_rawdata, i)
        command = ["bwa", "mem", reference_genome_fasta, fastq_file, "-t", str(threads) ]
        # redirect the output to the SAM file
        output_sam = os.path.join(os.getcwd(), "results_rna_seq",
                                  my_experiment, folder_results, i)  +".sam"
        
        with open(output_sam, "w") as output_file:
            subprocess.run(command, stdout=output_file, check=True)
        
def align_BWA_PE(fastq_files, path_rawdata, folder_results, reference_genome_fasta, end, threads):
    try:
      os.mkdir(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, folder_results))
    except OSError as error:
        print(error)
    experiments = []
    for file in fastq_files:
        parts = file.split(end)
        name = parts[0]
        experiments.append(name)
    
    experiments = set(experiments) 
    
    for i in experiments:  
        # Comamand for BWA
        input_fastq_r1 = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, path_rawdata, i+"_1.fastq.gz")
        input_fastq_r2 = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, path_rawdata, i+"_2.fastq.gz")
        command_BWA = ["bwa", "mem", "-t", "18", reference_genome_fna, input_fastq_r1, input_fastq_r2]
        # redirect the output to the SAM file
        output_sam = os.path.join(cur_dir, "results_rna_seq",
                                  my_experiment, folder_results, i)  +".sam"
        
        with open(output_sam, "w") as output_file:
            subprocess.run(command_BWA, stdout=output_file, check=True)


def generate_index_genome_bowtie(reference_fasta, reference_genome):
    bowtie_build_cmd = ["bowtie-build", reference_fasta, reference_genome]
    subprocess.run(bowtie_build_cmd, check=True)
    
def align_bowtie_PE(fastq_files, path_rawdata, folder_results, reference_genome, end, threads):
    try:
      os.mkdir(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, folder_results))
    except OSError as error:
        print(error)
    
    experiments = []
    for file in fastq_files:
        parts = file.split(end)
        name = parts[0]
        experiments.append(name)
    
    experiments = set(experiments) 
    
    for i in experiments:  
       # Comamand for 
       input_fastq_r1 = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, path_rawdata, i+"_1.fastq.gz")
       input_fastq_r2 = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, path_rawdata, i+"_2.fastq.gz")
       output_sam = os.path.join(os.getcwd(), "results_rna_seq",
                                 my_experiment, folder_results, i)  +".sam"
       bowtie_align_cmd = ["bowtie", "-t", "-p", str(threads), "-x", reference_genome, "-1", input_fastq_r1, 
                            "-2", input_fastq_r2, "-S", output_sam]
       subprocess.run(bowtie_align_cmd, check=True)


bam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "Bowtie2")
bam_files = [bam for bam in os.listdir(bam_folder) if (bam.endswith(".bam") &  bam.startswith("sorted"))]  

def HTSeq_count_SE(bam_files, bam_folder, folder_results, threads, feature_type, IDATTR, GFF_file):
    try:
      os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, folder_results))
    except OSError as error:
        print(error)

    for b in bam_files:  
        bam_file =  os.path.join(cur_dir, "results_rna_seq", my_experiment, bam_folder, b)
        output_htseq = os.path.join(cur_dir, "results_rna_seq",
                                  my_experiment, folder_results , b)  +".counts"
    
        #command_htseq = ["htseq-count", "-f", "bam", "-t", 
         #          "gene", "-i", "gene_id", "-n", "20", "-r", "pos", bam_file, reference_genome_gtf]
    
        #SE data
        command_htseq = ["htseq-count", "-f", "bam", "-t", feature_type, 
                         "-i", IDATTR, "-n", str(threads),
                         bam_file, GFF_file]
    
        # Ejecutar el comando y redirigir la salida al archivo de conteo
        with open(output_htseq, "w") as output:
            subprocess.run(command_htseq, stdout=output, check=True)

def HTSeq_count_PE(bam_files, bam_folder, folder_results, threads, feature_type, IDATTR, GFF_file):
    try:
        os.mkdir(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, folder_results))
    except OSError as error:
        print(error)
    
    for b in bam_files:  
        bam_file =  os.path.join(os.getcwd(), "results_rna_seq", my_experiment, bam_folder, b)
        output_htseq = os.path.join(os.getcwd(), "results_rna_seq",
                                 my_experiment, folder_results , b)  +".counts"
    
        command_htseq = ["htseq-count", "-f", "bam", "-t", feature_type, 
                         "-i", IDATTR, "-n", str(threads), 
                         "-r", "pos", bam_file, GFF_file]
    
    
    
        # Ejecutar el comando y redirigir la salida al archivo de conteo
        with open(output_htseq, "w") as output:
            subprocess.run(command_htseq, stdout=output, check=True)
            
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
ax.scatter(basemean_non_DE, log2_FC_non_DE, color='black', latype_databel = "non DE")
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


#%% #Extract fasta sequence of each replicon




from BCBio import GFF
from Bio import SeqIO

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

#%%  APERO data analysis


import os
import pandas as pd

def process_files(file1, file2, apero_folder):
    # Leer el primer archivo CSV y filtrar según los criterios dados
    df1 = pd.read_csv(os.path.join(apero_folder, file1))
    df1_coverage = df1[df1['freq'] >= 10]
    df1_internal = df1_coverage[df1_coverage['Class'] != 'Ai']
    
    # Leer el segundo archivo CSV y filtrar según los criterios dados
    df2 = pd.read_csv(os.path.join(apero_folder, file2))
    df2_coverage = df2[df2['freq'] >= 10]
    df2_internal = df2_coverage[df2_coverage['Class'] != 'Ai']
    
    # Los small RNAs son menores que 600 y además en este caso mi pipeline
    # para detectar sRNAs solo considera aquellos menores que esta distancia
    df1_lg = df1_internal[df1_internal['lg'] <= 600]
    df2_lg = df2_internal[df2_internal['lg'] <= 600]
    
    return df1_lg, df2_lg

def merge_rows(df_internal):
    df_sorted = df_internal.sort_values(by='Position').reset_index(drop=True)
    merged_rows = []
    start_index = 0
    current_index = 1
    max_lg_row = df_sorted.iloc[start_index]

    while current_index < len(df_sorted):
        position_diff = df_sorted.loc[current_index, 'Position'] - df_sorted.loc[current_index - 1, 'Position']

        if position_diff < 10:
            if df_sorted.loc[current_index, 'lg'] > max_lg_row['lg']:
                max_lg_row = df_sorted.iloc[current_index]
        else:
            merged_rows.append(max_lg_row)
            start_index = current_index
            max_lg_row = df_sorted.iloc[start_index]

        current_index += 1

    merged_rows.append(max_lg_row)
    result = pd.DataFrame(merged_rows)
    return result

def merge_results(results1, results2):
    merged_results = pd.DataFrame()

    for index, row1 in results1.iterrows():
        position1 = row1['Position']
        lg1 = row1['lg']

        matching_rows = results2[
            (results2['Position'] >= (position1 - 10)) &
            (results2['Position'] <= (position1 + 10))
        ]

        if not matching_rows.empty:
            max_lg_row = matching_rows.loc[matching_rows['lg'].idxmax()]

            if lg1 > max_lg_row['lg']:
                merged_results = pd.concat([merged_results, row1.to_frame().T])
            else:
                merged_results = pd.concat([merged_results, max_lg_row.to_frame().T])

    merged_results.reset_index(drop=True, inplace=True)
    return merged_results

def main(my_experiment, initial_string):
    cur_dir = os.getcwd()
    apero_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "APERO")

    replicates1 = [file for file in os.listdir(apero_folder) if (file.endswith("-1.csv") and file.startswith(initial_string))]
    replicates2 = [file for file in os.listdir(apero_folder) if (file.endswith("-2.csv") and file.startswith(initial_string))]

    replicates1 = sorted(set(replicates1))
    replicates2 = sorted(set(replicates2))

    processed_dataframes = []

    for file1, file2 in zip(replicates1, replicates2):
        try:
            df1_internal, df2_internal = process_files(file1, file2, apero_folder)
            processed_dataframes.append((df1_internal, df2_internal))
        except Exception as e:
            print(f"Error processing files {file1} and {file2}: {e}")

    dataframes_filtered_size = [element for element in processed_dataframes if all(len(df) > 0 for df in element)]
    
    merged_rows = []
    for i in dataframes_filtered_size:
        result_i = [merge_rows(df) for df in i]
        merged_rows.append(result_i)
    
    merged_replicates = []
    for i in merged_rows:
        merged_replicates.append(merge_results(i[0], i[1]))

    results_concatenated = pd.concat(merged_replicates, ignore_index=True)
    final_result = merge_rows(results_concatenated)
    return final_result

def apero_to_GFF(sRNAs_apero, seqid):
    gff = pd.DataFrame()
    gff['Column1'] = [seqid] * len(sRNAs_apero)
    gff['Column2'] = ['APERO'] * len(sRNAs_apero)
    sequential_numbers = range(1, len(sRNAs_apero) + 1)
    gff['Column3'] = ['predicted_sRNA_{}'.format(i) for i in sequential_numbers]
    gff['Column4'] = sRNAs_apero['Position']
    gff['Column5'] = sRNAs_apero['Position'] + sRNAs_apero['lg']
    gff['Column6'] = ['.'] * len(sRNAs_apero)
    gff['Column7'] = sRNAs_apero['str']
    gff['Column8'] = ['.'] * len(sRNAs_apero)
    gff['Column9'] = 'name=' + gff['Column3'] + '**' + sRNAs_apero['Position'].astype(str) + '-' + (sRNAs_apero['Position'] + sRNAs_apero['lg']).astype(str)
    gff['Column3'] = ['putative_sRNA'] * len(sRNAs_apero)
    return gff

def concat_gff(path_original, gff_apero):
    gff_original = pd.read_csv(path_original, comment='#', sep='\t', header=None)
    columnas_originals = gff_apero.columns
    gff_original.columns = columnas_originals
    gff_concat = pd.concat([gff_apero, gff_original])
    concatenated_sorted = gff_concat.sort_values(by=gff_concat.columns[3])
    return concatenated_sorted


my_experiment = "S-clav_RNA-seq_sRNas_OMEGA" 
sRNAs_PLASMID = main(my_experiment, "PLASMID")
sRNAs_PLASMID.reset_index(drop=True, inplace=True)
sRNAs_CHR = main(my_experiment, "CHR")
sRNAs_CHR.reset_index(drop=True, inplace=True)

gff_PLASMID = apero_to_GFF(sRNAs_PLASMID, "gi|99999999|ref|NC_999999.9|NZ_CP027859.1")
gff_CHR = apero_to_GFF(sRNAs_CHR, "gi|99999999|ref|NC_999999.9|NZ_CP027858.1")


df1_IGRs = sRNAs_PLASMID[sRNAs_PLASMID['Class'] == 'O']
df12_IGRs = sRNAs_CHR[sRNAs_CHR['Class'] == 'O']


frecuencia = sRNAs_CHR['Class'].value_counts()

path_CHR_original = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/NZ_CP027858.1.gff"
path_PLASMID_original = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/NZ_CP027859.1.gff"


#CHR
CHR_concatenated = concat_gff(path_CHR_original, gff_CHR)
CHR_concatenated['Column4'] = CHR_concatenated['Column4'].astype(int)
CHR_concatenated['Column5'] = CHR_concatenated['Column5'].astype(int)

#PLASMID
PLASMID_concatenated = concat_gff(path_PLASMID_original, gff_PLASMID)
PLASMID_concatenated['Column4'] = PLASMID_concatenated['Column4'].astype(int)
PLASMID_concatenated['Column5'] = PLASMID_concatenated['Column5'].astype(int)

#Save the results
file_gff_APERO_CHR = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "APERO", "CHR.gff")
with open(file_gff_APERO_CHR, 'w') as f:
    f.write('# produced by APERO\n')

# Save the DataFrame without indexes or column names
CHR_concatenated.to_csv(file_gff_APERO_CHR, sep='\t', index=False, header=False, mode='a')


file_gff_APERO_PLASMID = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "APERO", "PLASMID.gff")
with open(file_gff_APERO_PLASMID, 'w') as f:
    f.write('# produced by APERO\n')

# Save the DataFrame without indexes or column names
PLASMID_concatenated.to_csv(file_gff_APERO_PLASMID, sep='\t', index=False, header=False, mode='a')


#%% Counts reads alignment with HTseq-count once the 
#sRNAshas been predicted 
bam_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_CHR")
bam_files = [bam for bam in os.listdir(bam_folder) if (bam.endswith(".bam") &  bam.startswith("sorted"))]  
file_gff_APERO_CHR = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "APERO", "CHR.gff")

HTSeq_count_PE(bam_files, bam_folder, "HTseq_APERO_CHR", 18, "putative_sRNA", "name", file_gff_APERO_CHR)


#Correr para el plasmido
bam_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_PLASMID")
bam_files = [bam for bam in os.listdir(bam_folder) if (bam.endswith(".bam") &  bam.startswith("sorted"))]  
file_gff_APERO_PLASMID = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "APERO", "PLASMID.gff")


HTSeq_count_PE(bam_files, bam_folder, "HTseq_APERO_PLASMID", 18, "putative_sRNA", "name", file_gff_APERO_PLASMID)

#%% use deseq2 to determine if the genes are differentially expressed

def DESeq2(results_folder, HTSEQ_folder, desired_I, MD_file):
    try:
      os.mkdir(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, results_folder))
    except OSError as error:
        print(error)
    DESeq2_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, results_folder)


    #Convert the data from HTseq to the format required by DESeq2
    HTseq_folder = HTSEQ_folder
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

    desired_index = desired_I
    HT_concat= HT_concat.reindex(index=desired_index)

    #save the results
    HT_file = os.path.join(DESeq2_folder, "HTseq.txt")
    HT_concat.to_csv(HT_file, sep = "\t", index=True)

    #Open the metadata file
    metadata_file = MD_file
    metadata = pd.read_table(metadata_file, sep = "\t")
    metadata.set_index('Sample', inplace=True)


    #Filter out genes that have less than 5 read counts in total
    genes_to_keep = HT_concat.columns[HT_concat.sum(axis=0) >= 10]
    HT_concat = HT_concat[genes_to_keep]



    # DeseqDataSet fits dispersion and log-fold change (LFC) parameters from
    #the data, and stores them
    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=HT_concat,
        metadata=metadata,
        design_factors="Condition",
        refit_cooks=True,
        inference=inference,)
    return dds

DESeq2_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "DESeq2")
results_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "DESeq2_APERO_CHR")
HTseq_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "HTseq_APERO_CHR")
desired_index = ['SC24H-1', 'SC24H-2', 'SC48H-1', 'SC48H-2', 'SC72H-1', 'SC72H-2', 'SC96H-1', 'SC96H-2']
metadata_file = DESeq2_folder+"/metadata.txt"


deseq2_CHR = DESeq2(results_folder, HTseq_folder, desired_index, metadata_file)


#Once a DeseqDataSet was initialized, we may run the deseq2() method to 
#fit dispersions and LFCs.
deseq2_CHR.deseq2()
    
#save the results
with open(os.path.join(DESeq2_folder, "dds_detailed_pipe.pkl"), "wb") as f:
            pkl.dump(dds, f)

#Now that dispersions and LFCs were fitted, we may proceed with
#statistical tests to compute p-values and adjusted p-values for 
#differential expresion. 
inference = DefaultInference(n_cpus=8)
stat_res = DeseqStats(deseq2_CHR, inference=inference)
    
#PyDESeq2 computes p-values using Wald tests. 
#This can be done using the summary() method
stat_res.summary()
   
summary_DESeq2= stat_res.results_df
    
#For visualization or post-processing purposes, it might be suitable
#to perform LFC shrinkage. 
lfc_shrink = stat_res.lfc_shrink(coeff="Condition_Late-exponential_vs_Early-exponential")
#%% plots
## PCA
sc.tl.pca(deseq2_CHR)
sc.pl.pca(deseq2_CHR, color = 'Condition', size = 200)

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