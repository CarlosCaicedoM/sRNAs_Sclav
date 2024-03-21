#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 08:11:39 2024

@author: Carlos_Caicedo-Montoya
"""

import subprocess
import os
import gzip
import re

#Create file to download your data
my_experiment = "S-clavuligerus_RNA-Seq-data_SRA"
my_file = "Sclav_RNA-seq_data.txt"

reference_genome_fna = "/home/usuario/Documentos/Paula_PhD/raw_data/S-clavuligerus_ATCC27064_annotations/GCF_005519465.1_ASM551946v1_genomic.fna"
reference_genome_gtf = "/home/usuario/Documentos/Paula_PhD/raw_data/S-clavuligerus_ATCC27064_annotations/GCF_005519465.1_ASM551946v1_genomic.gtf"
reference_genome = "SCLAV"

reference_chr_fasta = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/NZ_CP027858.1/NZ_CP027858.1.fna"
reference_plasmid_fasta = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/NZ_CP027859.1/NZ_CP027859.1.fna"
reference_chr = "SCLAV_CHR"
reference_plasmid = "SCLAV_PLASMID"


cur_dir = os.getcwd()
path_rawdata = os.path.join(os.getcwd(), "raw_data", my_experiment) 


try: 
    os.mkdir(path_rawdata) 
except OSError as error: 
    print(error) 

my_data = os.path.join(os.getcwd(), my_file)
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


# List to store the names of the fastq files
fastq_files = []

# Read all files in the downloads directory
fastq_files = [file for file in os.listdir(path_rawdata) if file.endswith(".fastq") ]

# Compress the fastq files and delete the originals
for filename in fastq_files:
    # Path to the original fastq file
    original_path = os.path.join(path_rawdata, filename)  
    # Path to the compressed file
    compressed_path = os.path.join(path_rawdata, filename + ".gz")
    
    # Open the original fastq file in binary read mode
    with open(original_path, 'rb') as f_in:
        # Open the compressed file in binary write mode
        with gzip.open(compressed_path, 'wb') as f_out:
            # Copy data from the original fastq file to the compressed file
            f_out.writelines(f_in)

    # Delete the original fastq file
    os.remove(original_path)

import os
import gzip
from concurrent.futures import ThreadPoolExecutor

def compress_fastq(filename):
    # Path to the original fastq file
    original_path = os.path.join(path_rawdata, filename)  
    # Path to the compressed file
    compressed_path = os.path.join(path_rawdata, filename + ".gz")
    
    # Open the original fastq file in binary read mode
    with open(original_path, 'rb') as f_in:
        # Open the compressed file in binary write mode
        with gzip.open(compressed_path, 'wb') as f_out:
            # Copy data from the original fastq file to the compressed file
            f_out.writelines(f_in)


max_workers = 18  # Por ejemplo, 4 hilos

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    executor.map(compress_fastq, fastq_files)


#%% Initial quality control of the data with fastqc
try:
    os.mkdir("results_rna_seq")
except OSError as error:
    print(error)  

try:
    os.mkdir(os.path.join("results_rna_seq", my_experiment))
except OSError as error:
    print(error)
fastqc_outdir=os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "fastqc")
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

MultiQC(fastqc_outdir, "multiqc")

#%% Align sequences to the reference genome  Bowtie2

try:
  os.mkdir(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2"))
except OSError as error:
    print(error)



# generate the index genome file

# Comando Bowtie2 para crear la base de datos de referencia
bowtie2_build_cmd = ["bowtie2-build", reference_genome_fna, reference_genome]
subprocess.run(bowtie2_build_cmd, check=True)


#SE data
fastq_files_compressed = [file for file in os.listdir(path_rawdata) if file.endswith(".gz") ]

fastq_files_compressed.remove('SRR12414212_1.fastq.gz')
fastq_files_compressed.remove('SRR12414212_2.fastq.gz')

for i in fastq_files_compressed:  
    fastq_file = os.path.join(path_rawdata, i)
    output_sam = os.path.join(os.getcwd(), "results_rna_seq",
                              my_experiment, "Bowtie2", i)  +".sam"
    bowtie2_align_cmd = ["bowtie2", "-p", "18",
                         "-x",  reference_genome,
                         "-U", fastq_file,
                         "-S",  output_sam]
    subprocess.run(bowtie2_align_cmd, check=True)

#PE_DATA
fastq_files_PE = ['SRR12414212_1.fastq.gz', 'SRR12414212_2.fastq.gz']
input_fastq_r1 = os.path.join(path_rawdata,  fastq_files_PE[0])
input_fastq_r2 = os.path.join(path_rawdata,  fastq_files_PE[1])

output_sam = os.path.join(os.getcwd(), "results_rna_seq",
                              my_experiment, "Bowtie2", "SRR12414212")  +".sam"
bowtie2_align_cmd = ["bowtie2", "-p", "19", "-x", reference_genome, 
                     "-1", input_fastq_r1 , 
                     "-2", input_fastq_r2 , "-S", output_sam]
subprocess.run(bowtie2_align_cmd, check=True)

#%% Convert SAM files to BAM 
import pysam


alignment_sam_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2")
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
sorted_bam(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2"))

#%%  Delete sam files


def delete_files_aln(folder, suffix):
    files = os.listdir(folder)
    for file in files:
        path_file = os.path.join(folder, file)
        if file.endswith(suffix) or (file.endswith(".bam") and not file.startswith("sorted")):
            os.remove(path_file)


delete_files_aln(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2"), ".sam")


#%% Align sequences to the chromosome and plasmid individually with with  Bowtie2

# generate the index genome file
def generate_index_genome(reference_fasta, reference_genome):
    bowtie2_build_cmd = ["bowtie2-build", reference_fasta, reference_genome]
    subprocess.run(bowtie2_build_cmd, check=True)

generate_index_genome(reference_chr_fasta, reference_chr)
generate_index_genome(reference_plasmid_fasta, reference_plasmid)



#SE data
pattern = re.compile(r"^(?!.*(?:_1|_2)\.(?:fastq\.gz|gz|fastq))")
fastq_files_compressed_SE = [file for file in os.listdir(path_rawdata) if pattern.match(file)]



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

align_bowtie2_SE(fastq_files_compressed_SE, path_rawdata, "Bowtie2_CHR", reference_chr, 18)
align_bowtie2_SE(fastq_files_compressed_SE, path_rawdata, "Bowtie2_PLASMID", reference_plasmid, 18)

#PE_DATA
pattern = re.compile(r".*_[12]\.(?:fastq\.gz|gz|fastq)$")
fastq_files_compressed_PE = [file for file in os.listdir(path_rawdata) if pattern.match(file)]

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

align_bowtie2_PE(fastq_files_compressed_PE, path_rawdata, "Bowtie2_CHR", reference_chr, "_", 18)
align_bowtie2_PE(fastq_files_compressed_PE, path_rawdata, "Bowtie2_PLASMID", reference_plasmid, "_", 18)

    
#%% Convert SAM files to BAM files 
  
# Alignment for PE CHR
alignment_sam_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_CHR")
alignment_sam = [sam for sam in os.listdir(alignment_sam_folder) if sam.endswith(".sam") ]  
for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    output_bam = os.path.join(alignment_sam_folder, os.path.splitext(alignment)[0]+".bam")
    convert_sam_to_bam(input_sam, output_bam)

  
# Alignment for PE plasmid
alignment_sam_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_PLASMID")
alignment_sam = [sam for sam in os.listdir(alignment_sam_folder) if sam.endswith(".sam") ]  
for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    output_bam = os.path.join(alignment_sam_folder, os.path.splitext(alignment)[0]+".bam")
    convert_sam_to_bam(input_sam, output_bam)
#%% Convert unsorted to sorted BAM files

sorted_bam(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_CHR"))
sorted_bam(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_PLASMID"))

#%%  Delete sam files

delete_files_aln(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_CHR"), ".sam")
delete_files_aln(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_PLASMID"), ".sam")

#%% generate .bai files

#access a BAM index file
def bam_index_file(alignment_sam_folder, alignment_bam_sorted):
    for i in alignment_bam_sorted:
        bam_file = os.path.join(alignment_sam_folder, i)
        command_index = ["samtools", "index", bam_file]    
        subprocess.run(command_index)

#CHR
alignment_sam_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_CHR")
alignment_bam_sorted = [bam for bam in os.listdir(alignment_sam_folder) if bam.startswith("sorted")]  
bam_index_file(alignment_sam_folder, alignment_bam_sorted)


#PLASMID
alignment_sam_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_PLASMID")
alignment_bam_sorted = [bam for bam in os.listdir(alignment_sam_folder) if bam.startswith("sorted")]  
bam_index_file(alignment_sam_folder, alignment_bam_sorted)


#%% Align data Luisa


os.chdir("/home/usuario/Documentos/Luisa_PhD")
my_experiment = "S-clav_mutants_RNA-seq"
fastq_filtered_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "fastq_rRNA_removed")
fastq_filtered = [fastq for fastq in os.listdir(fastq_filtered_folder) if fastq.endswith(".gz") ]  
fastq_filtered.remove("PACP72H-1_1.fastq.gz")
fastq_filtered.remove("PACP72H-1_2.fastq.gz")
fastq_filtered.remove("PACP72H-2_1.fastq.gz")
fastq_filtered.remove("PACP72H-2_2.fastq.gz")


generate_index_genome(reference_chr_fasta, reference_chr)
generate_index_genome(reference_plasmid_fasta, reference_plasmid)

align_bowtie2_PE(fastq_filtered, fastq_filtered_folder, "Bowtie2_CHR", reference_chr, "_", 18)
align_bowtie2_PE(fastq_filtered, fastq_filtered_folder, "Bowtie2_PLASMID", reference_plasmid, "_", 18)

# Alignment for PE CHR
alignment_sam_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_CHR")
alignment_sam = [sam for sam in os.listdir(alignment_sam_folder) if sam.endswith(".sam") ]  
for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    output_bam = os.path.join(alignment_sam_folder, os.path.splitext(alignment)[0]+".bam")
    convert_sam_to_bam(input_sam, output_bam)

  
# Alignment for PE plasmid
alignment_sam_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_PLASMID")
alignment_sam = [sam for sam in os.listdir(alignment_sam_folder) if sam.endswith(".sam") ]  
for alignment in alignment_sam: 
    input_sam = os.path.join(alignment_sam_folder, alignment)
    output_bam = os.path.join(alignment_sam_folder, os.path.splitext(alignment)[0]+".bam")
    convert_sam_to_bam(input_sam, output_bam)
#Convert unsorted to sorted BAM files

sorted_bam(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_CHR"))
sorted_bam(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_PLASMID"))

#Delete sam files

delete_files_aln(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_CHR"), ".sam")
delete_files_aln(os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_PLASMID"), ".sam")

#CHR
alignment_sam_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_CHR")
alignment_bam_sorted = [bam for bam in os.listdir(alignment_sam_folder) if bam.startswith("sorted")]  
bam_index_file(alignment_sam_folder, alignment_bam_sorted)


#PLASMID
alignment_sam_folder = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "Bowtie2_PLASMID")
alignment_bam_sorted = [bam for bam in os.listdir(alignment_sam_folder) if bam.startswith("sorted")]  
bam_index_file(alignment_sam_folder, alignment_bam_sorted)
