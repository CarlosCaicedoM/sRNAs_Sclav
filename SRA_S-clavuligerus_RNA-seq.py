#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 08:11:39 2024

@author: Carlos_Caicedo-Montoya
"""

import subprocess
import os
import gzip

#Create file to download your data
my_experiment = "S-clavuligerus_RNA-Seq-data_SRA"
my_file = "Sclav_RNA-seq_data.txt"

cur_dir = os.getcwd()
path = os.path.join(cur_dir, "raw_data", my_experiment) 


try: 
    os.mkdir(path) 
except OSError as error: 
    print(error) 

my_data = os.path.join(cur_dir, my_file)
#If you want to remove all whitespace characters (newlines and spaces) from the end of each line
with open(my_data) as f:
    lines = [line.rstrip() for line in f]
Accessions  = lines[0:-1]

# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in Accessions:
    print ("Currently downloading: " + sra_id)
    fasterq = "fasterq-dump " + sra_id + " -O " + path
    print ("The command used was: " + fasterq)
    subprocess.call(fasterq, shell=True)


# List to store the names of the fastq files
fastq_files = []

# Read all files in the downloads directory
fastq_files = [file for file in os.listdir(path) if file.endswith(".fastq") ]

# Compress the fastq files and delete the originals
for filename in fastq_files:
    # Path to the original fastq file
    original_path = os.path.join(path, filename)  
    # Path to the compressed file
    compressed_path = os.path.join(path, filename + ".gz")
    
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
    original_path = os.path.join(path, filename)  
    # Path to the compressed file
    compressed_path = os.path.join(path, filename + ".gz")
    
    # Open the original fastq file in binary read mode
    with open(original_path, 'rb') as f_in:
        # Open the compressed file in binary write mode
        with gzip.open(compressed_path, 'wb') as f_out:
            # Copy data from the original fastq file to the compressed file
            f_out.writelines(f_in)

# Suponiendo que `fastq_files` es una lista de los nombres de los archivos FASTQ
# y `path` es la ruta donde se encuentran los archivos.
# Puedes ajustar `max_workers` según la cantidad de hilos que desees utilizar.
max_workers = 18  # Por ejemplo, 4 hilos

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    executor.map(compress_fastq, fastq_files)


