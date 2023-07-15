#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 13:22:39 2022

@author: Carlos Caicedo-Montoya
"""


#conda install -c bioconda fastqc
#sudo apt install sra-toolkit
#conda install -c bioconda -c conda-forge multiqc
#conda create -n cutadaptenv cutadapt
#sudo apt install cutadapt
# conda install -c bioconda fastp
#sudo apt install bwa
#sudo apt install bowtie2

#pip install HTSeq
#pip install pydeseq2
#conda install -c conda-forge scanpy python-igraph leidenalg
#pip install gseapy
#pip install goatools

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
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import scanpy as sc
import gseapy as gp
from gseapy.plot import gseaplot
import numpy as np
import matplotlib.pyplot as plt


#%% get the data neccesary to run this pipeline
#Create file to download your data
#my_experiment = "S-enterica_APERO"
my_experiment = "S-clavuligerus_Hwang-2019"
#my_file = "/SraAccList_S-enterica_APERO.txt"
my_file = "SraAccList_S_clavuligerus_Hwang_2019.txt"
#reference_genome = "/media/usuario/CACM/Documentos/sRNAs_Sclav/raw_data/S-enterica_APERO/GCF_000210855.2_ASM21085v2_genomic.fna"
reference_genome = "/home/usuario/Documentos/sRNAS_Sclav/raw_data/S-clavuligerus_ATCC27064_annotations/GCF_005519465.1_ASM551946v1_genomic.fna"

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
	
fastq_files = [file for file in os.listdir(path) if file.endswith(".fastq") ]
fastq_files_path = []
for i in fastq_files:
    fastq_files_path.append(os.path.join(path, i))

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

fastq_command =[]  
for i in fastq_files:
    fastq_command.append(os.path.join(path, i))
fastq_command.insert(0, "fastqc")
fastq_command.append("-o")
fastq_command.append(fastqc_outdir)
fastq_command.append("-t")
fastq_command.append("20")   #change the number of threads accordingly
#--noextract
subprocess.call(fastq_command)

#Summarize all the fastqc results with multiQC
try:
    os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "multiqc"))
except OSError as error:
    print(error)
 
multiQC_command = ["multiqc", fastqc_outdir]
multiQC_command.append("-n")
multiQC_command.append(my_experiment)
multiQC_command.append("-o")
multiQC_command.append(os.path.join(cur_dir, "results_rna_seq", my_experiment, "multiqc"))
multiQC_command.append("-f")
multiQC_command.append("--profile-runtime") 
subprocess.call(multiQC_command)

#%% Filter the adpaters with cutadapt
"""
try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "cutadapt"))
except OSError as error:
    print(error) 

adapter_5prime = "GTTCAGAGTTCTACAGTCCGACGATC"
adapter_3prime = "TGGAATTTCTCGGGTGCCAAGG"


#Cutadapt command for adapter trimming at both ends, 
#quality triming at both ends and filtering to only have reads with length >18nt


for i,j  in enumerate(fastq_files):
    cutadapt_command=["cutadapt", "-a", adapter_3prime, "-g",
                  adapter_5prime, "-q", "15,10", "-m", "18", 
                  "--cores=0", "-o",  
                  os.path.join(cur_dir, "results_rna_seq", my_experiment, "cutadapt", "trimmed_"+j), fastq_files_path[i]]

    subprocess.call(cutadapt_command)
"""   

#%% Fastqc after cutadapt


fastq_command_cutadapt =[]  
for i, j  in enumerate(fastq_files):
    fastq_command_cutadapt.append(os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+j))
fastqc_trimmed=os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastqc_cutadapt")
try:
    os.mkdir(fastqc_trimmed)
except OSError as error:
    print(error) 

   
fastq_command_cutadapt.insert(0, "fastqc")
fastq_command_cutadapt.append("-o")
fastq_command_cutadapt.append(fastqc_trimmed)
fastq_command_cutadapt.append("-t")
fastq_command_cutadapt.append("20")   #change the number of threads accordingly
#--noextract
subprocess.call(fastq_command_cutadapt)

#Summarize all the fastqc results with multiQC
try:
    os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "multiqc_cutadapt"))
except OSError as error:
    print(error)
 
multiQC_command = ["multiqc", fastqc_trimmed]
multiQC_command.append("-n")
multiQC_command.append(my_experiment)
multiQC_command.append("-o")
multiQC_command.append(os.path.join(cur_dir, "results_rna_seq", my_experiment, "multiqc_cutadapt"))
multiQC_command.append("-f")
multiQC_command.append("--profile-runtime") 
subprocess.call(multiQC_command)

#%% Quality control and filtering by fastp
try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp"))
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

for i,j  in enumerate(Accessions):
    fastp_command=["fastp", 
                   "-i", path+"/"+j+"_1.fastq",
                   "-I",  path+"/"+j+"_2.fastq",               
                   "-o", os.path.join(cur_dir, "results_rna_seq",
                                my_experiment, "fastp", "trimmed_"+j+"_1.fastq"), 
                   "-O", os.path.join(cur_dir, "results_rna_seq",
                                my_experiment, "fastp", "trimmed_"+j+"_2.fastq"),
                   "-j", j+".json", "-h", j+".html", "-w", "16", "--cut_right"]

    subprocess.call(fastp_command)

#Move the reports to the fastp folder
reports = [rep for rep in os.listdir(os.curdir) if (rep.endswith(".html") or rep.endswith(".json"))]  
for file_name in reports:
    shutil.move(file_name, os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp"))




#adapter trimming is enabled by default
#quality filtering is enabled by default. 
#length filtering is enabled by default.
#-w number of threads 
#-detect_adapter_for_pe        
#by default, the auto-detection for adapter is for SE data input only,
#turn on this option to enable it for PE data

"""
-r, --cut_right move a sliding window from front to tail, 
if meet one window with mean quality < threshold, drop the bases in the 
window and the right part, and then stop. Use cut_right_window_size to set 
the widnow size, and cut_right_mean_quality to set the mean quality threshold. 
This is similar as the Trimmomatic SLIDINGWINDOW method
"""


#%% Fastqc after trimming

fastq_command_fastp =[]  
for i, j  in enumerate(fastq_files):
    fastq_command_fastp.append(os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", "trimmed_"+j))
fastqc_trimmed=os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastqc_trimmed")
try:
    os.mkdir(fastqc_trimmed)
except OSError as error:
    print(error) 

   
fastq_command_fastp.insert(0, "fastqc")
fastq_command_fastp.append("-o")
fastq_command_fastp.append(fastqc_trimmed)
fastq_command_fastp.append("-t")
fastq_command_fastp.append("20")   #change the number of threads accordingly
#--noextract
subprocess.call(fastq_command_fastp)

#Summarize all the fastqc results with multiQC
try:
    os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "multiqc_fastp"))
except OSError as error:
    print(error)
 
multiQC_command = ["multiqc", fastqc_trimmed]
multiQC_command.append("-n")
multiQC_command.append(my_experiment)
multiQC_command.append("-o")
multiQC_command.append(os.path.join(cur_dir, "results_rna_seq", my_experiment, "multiqc_fastp"))
multiQC_command.append("-f")
multiQC_command.append("--profile-runtime") 
subprocess.call(multiQC_command)

#%% Align sequences to the reference genome


try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "BWA"))
except OSError as error:
    print(error)


fastq_filtered_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp")
fastq_filtered = [fastq for fastq in os.listdir(fastq_filtered_folder) if fastq.endswith(".fastq") ]  

# generate the index genome file
command = ['bwa', 'index', reference_genome]
subprocess.run(command, check=True)


#SE data
for i in fastq_filtered:  
    # Comamand for BWA
    fastq_file = os.path.join(cur_dir, "results_rna_seq", my_experiment, "fastp", i)
    command = ["bwa", "mem", reference_genome, fastq_file]
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




#%% count features with htseq-count

try:
  os.mkdir(os.path.join(cur_dir, "results_rna_seq", my_experiment, "HTseq"))
except OSError as error:
    print(error)


reference_genome_gtf= "/media/usuario/CACM/Documentos/sRNAs_Sclav/raw_data/annotations_Streptococcus/GCF_000007465.2_ASM746v2_genomic.gtf"
sam_folder = os.path.join(cur_dir, "results_rna_seq", my_experiment, "BWA")
sam_files = [sam for sam in os.listdir(sam_folder) if sam.endswith(".sam") ]  


for s in sam_files:  
    # Definir los archivos de entrada y salida
    sam_file =  os.path.join(cur_dir, "results_rna_seq", my_experiment, "BWA", s)
    output_htseq = os.path.join(cur_dir, "results_rna_seq",
                              my_experiment, "HTseq", s)  +".counts"

    # Comando para HTSeq-count
    command_htseq = ["htseq-count", "-f", "sam", "-t", 
               "gene", "-i", "gene_id", sam_file, reference_genome_gtf]

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
    match = re.search(r"SRR\d+", HT)
    if match:
        resultado = match.group(0)
    df = df.rename(columns={1: resultado})
    df.set_index('Gene_ID', inplace=True)
    dfs_HTseq.append(df)

HT_concat = pd.concat(dfs_HTseq, axis = 1)

rows_eliminate = df.index[-5:]  # Get the index of the last 5 rows
HT_concat.drop(rows_eliminate, inplace=True)

HT_concat = HT_concat.transpose()
HT_file = os.path.join(DESeq2_folder, "HTseq.txt")
HT_concat.to_csv(HT_file, sep = "\t", index=True)
metadata_file = DESeq2_folder+"/metadata.txt"
metadata = pd.read_table(metadata_file, sep = "\t")
metadata.set_index('Sample', inplace=True)


#Filter out genes that have less than 10 read counts in total
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

#save the resukts
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
lfc_shrink = stat_res.lfc_shrink(coeff="Condition_Treatment_vs_Control")

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
FC_thr = 1.8
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


#%% Gene enrichment analysis

#%% Explore GSEApy

"""https://gseapy.readthedocs.io/en/latest/introduction.html"""
ranking = summary_DESeq2[['Gene_ID', 'stat']].dropna().sort_values('stat', ascending = False)




#%% como saber si los datos son stranded o no 



#%%Simple quality filtering for FASTQ files
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO

for sra_id in Accessions:
    fastq_file = sra_id+".fastq"
    count = 0
    #Count sequences
    for rec in SeqIO.parse(os.path.join(path, fastq_file), "fastq"):
        count += 1
    print("%i reads" % count)

    #Filter by quality 
    good_reads = (
        rec
    for rec in SeqIO.parse(os.path.join(path, fastq_file), "fastq")
        if min(rec.letter_annotations["phred_quality"]) >= 20
            )
    count = SeqIO.write(good_reads, "good_quality"+fastq_file, "fastq")
    print("Saved %i reads" % count)


#Trimming off primer sequences
count = 0
total_len = 0
with open(os.path.join(path, fastq_file)) as in_handle:
    for title, seq, qual in FastqGeneralIterator(in_handle):
        count += 1
        total_len += len(seq)
print("%i records with total sequence length %i" % (count, total_len))

#Sequences with the primers
primer_reads = (
    rec
    for rec in SeqIO.parse(os.path.join(path, fastq_file), "fastq")
    if rec.seq.startswith("GATGACGGTGT")
)
count = SeqIO.write(primer_reads, "with_primer.fastq", "fastq")
print("Saved %i reads" % count)


#Sequences with the primers but trimmed
trimmed_primer_reads = (rec[11:]
    for rec in SeqIO.parse(os.path.join(path, fastq_file), "fastq")
    if rec.seq.startswith("GATGACGGTGT"))
count = SeqIO.write(trimmed_primer_reads, "with_primer_trimmed.fastq", "fastq")
print("Saved %i reads" % count)



def trim_primer(record, primer):
    if record.seq.startswith(primer):
        return record[len(primer) :]
    else:
        return record


trimmed_reads = (trim_primer(record, "GATGACGGTGT")
    for record in SeqIO.parse(os.path.join(path, fastq_file), "fastq"))
count = SeqIO.write(trimmed_reads, "trimmed.fastq", "fastq")
print("Saved %i reads" % count)


def trim_primers(records, primer):
    """Removes perfect primer sequences at start of reads.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    len_primer = len(primer)  # cache this for later
    for record in records:
        if record.seq.startswith(primer):
            yield record[len_primer:]
        else:
            yield record


original_reads = SeqIO.parse(os.path.join(path, fastq_file), "fastq")
trimmed_reads = trim_primers(original_reads, "GATGACGGTGT")
count = SeqIO.write(trimmed_reads, "trimmed.fastq", "fastq")
print("Saved %i reads" % count)


#Trimming off adaptor sequences
#This time however, we will look for the sequence anywhere in the reads, 
#not just at the very beginning
def trim_adaptors(records, adaptor):
    """Trims perfect adaptor sequences.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    len_adaptor = len(adaptor)  # cache this for later
    for record in records:
        index = record.seq.find(adaptor)
        if index == -1:
            # adaptor not found, so won't trim
            yield record
        else:
            # trim off the adaptor
            yield record[index + len_adaptor :]
original_reads = SeqIO.parse(os.path.join(path, fastq_file), "fastq")
trimmed_reads = trim_adaptors(original_reads, "GATGACGGTGT")
count = SeqIO.write(trimmed_reads, "trimmed.fastq", "fastq")
print("Saved %i reads" % count)

#let’s add a minimum length requirement as well:
def trim_adaptors(records, adaptor, min_len):
    """Trims perfect adaptor sequences, checks read length.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    len_adaptor = len(adaptor)  # cache this for later
    for record in records:
        len_record = len(record)  # cache this for later
        if len(record) < min_len:
            # Too short to keep
            continue
        index = record.seq.find(adaptor)
        if index == -1:
            # adaptor not found, so won't trim
            yield record
        elif len_record - index - len_adaptor >= min_len:
            # after trimming this will still be long enough
            yield record[index + len_adaptor :]

original_reads = SeqIO.parse("SRR020192.fastq", "fastq")
trimmed_reads = trim_adaptors(original_reads, "GATGACGGTGT", 100)
count = SeqIO.write(trimmed_reads, "trimmed.fastq", "fastq")
print("Saved %i reads" % count)

