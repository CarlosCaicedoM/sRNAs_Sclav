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
import multiqc
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
path_rawdata = os.path.join(cur_dir, "RNA-seq_OMEGA")


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

fastQC_command =[]  
for i in fastq_files:
    fastQC_command.append(os.path.join(path_rawdata, i))
fastQC_command.insert(0, "fastqc")
fastQC_command.append("-o")
fastQC_command.append(fastqc_outdir)
fastQC_command.append("-t")
fastQC_command.append("20")   #change the number of threads accordingly
#--noextract
subprocess.call(fastQC_command)


