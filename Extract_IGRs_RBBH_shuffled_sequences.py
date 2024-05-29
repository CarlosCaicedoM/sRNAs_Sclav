#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 25 13:23:11 2024

@author: usuario
"""


#install all dependencies


#conda install -c conda-forge biopython 
#conda install -c bioconda blast
#conda install -c bioconda transtermhp 

#%%  Functions
from IPython import get_ipython
get_ipython().magic('reset -sf')
#run_line_magic(magic_name, parameter_s).get_ipython().magic('reset -sf')

import time
from random import randint
from Bio import SeqIO, SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
import os
from os import remove
import subprocess
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
import pandas as pd
from Bio import SearchIO
from Bio import AlignIO
import shutil
import matplotlib.pyplot as plt
import numpy as np
import collections
import re
from Bio.SeqUtils import gc_fraction
from collections import Counter
from scipy.stats import chi2_contingency
import random
from pandas.errors import EmptyDataError

#from Bio.SeqUtils import GC

### Define function to extract IGRs
def get_intergenic_regions(handle, intergenic_length, max_intergenic_length, output_dir = "results_sRNA/IGRs"):   

    intergenic_records = []
    
    for seq_record in SeqIO.parse(handle, "genbank"):
        gene_locations = []
        for feature in seq_record.features:
            if feature.type == 'gene':   #'gene'; 'CDS'
                my_start = feature.location.start
                my_end = feature.location.end
                my_strand = feature.strand
                gene_ID = str(feature.qualifiers['locus_tag'][0])
                gene_locations.append((my_start, my_end, my_strand, gene_ID))

        for i,location in enumerate(gene_locations[1:]):    
            # Compare current start position to previous end position
            last_end = gene_locations[i][1]
            this_start = location[0]
            strand = location[2]
            last_ID = gene_locations[i][3]
            current_ID = location[3]
            last_strand = gene_locations[i][2]
            if this_start > last_end and this_start - last_end >= intergenic_length and this_start - last_end <= max_intergenic_length:
                intergene_seq = seq_record.seq[last_end:this_start]               
            else:
                continue
                        
            if strand ==1 and last_strand==-1:
                orientation = "DP"
            if strand ==-1 and last_strand ==1:
                orientation = "DT"
            if strand ==1 and last_strand ==1:
                orientation = "CO_F"
            if strand ==-1 and last_strand==-1:
                orientation = "CO_R"
                        
            intergenic_records.append(
                              SeqRecord(intergene_seq,
                                        id="%s++%s-%s++%d-%d++%s" % (seq_record.id, last_ID, 
                                                                         current_ID,
                                                                         last_end,
                                                                         this_start,
                                                                         orientation), 
                                        description=""))
                               
    base_name=seq_record.dbxrefs[2]
    base_name_list=base_name.split(":")
    base_name_formated = base_name_list[1]
    outpath = output_dir + "/" + os.path.splitext(os.path.basename(base_name_formated))[0] + "_IGRs.fasta"
    SeqIO.write(intergenic_records, open(outpath,"w"), "fasta")            
 

def trim_sequences(my_file, left, right, output_dir = "results_sRNA/IGRs"):
    intergenic_records_reduced= []

    #left and right must be less than len(seq)/2
    #Right must be greater than zero
    
    for seq_record in SeqIO.parse(my_file, "fasta"):
        if right > 0:
            intergene_seq_reduced = seq_record.seq[left:-right]
        else:
            intergene_seq_reduced = seq_record.seq[left:]
        
        description = seq_record.description.split("++")
        
        location = description[2].split("-")
        
        start = int(location[0])+left
        end = int(location[1])-right

        intergenic_records_reduced.append(SeqRecord(intergene_seq_reduced,
                                        id="%s++%s++%d-%d++%s" % (description[0], 
                                                                  description[1],
                                                                  start, 
                                                                  end, 
                                                                  description[3]), 
                                        description=""))
                                                                                                                                    
    base_name_formated = my_file
    outpath = output_dir + "/" + os.path.splitext(os.path.basename(base_name_formated))[0] + "_reduced.fasta"
    SeqIO.write(intergenic_records_reduced, open(outpath,"w"), "fasta")   

def make_blast_db(output_dir = "results_sRNA/IGRs"):
    print("Making Blast Databases")
    my_results_reduced = os.listdir(output_dir)
    for IGRs_reduced in my_results_reduced:
        IGRs_database =  os.path.join(output_dir, IGRs_reduced)
        cline = NcbimakeblastdbCommandline(dbtype="nucl", 
                                           input_file=IGRs_database, 
                                           out = IGRs_database)
        stdout, stderr = cline()


#Make blast database from IGRs sequence to compare
def my_genome(my_reference, output_dir = "results_sRNA/IGRs"):
    my_new_reference  = os.path.join(my_input_dir, my_reference)
    for seq_record in SeqIO.parse(my_new_reference, "genbank"):
        base_name=seq_record.dbxrefs[2]
        base_name_list=base_name.split(":")
        base_name_formated = base_name_list[1]
        outpath = output_dir + "/" + os.path.splitext(os.path.basename(base_name_formated))[0] + "_IGRs_reduced.fasta"
        my_reference_IGRs = outpath
    return my_reference_IGRs


def my_subjects(my_reference, output_dir = "results_sRNA/IGRs"):
    my_reference_IGRs = my_genome(my_reference)
    my_subjects = []
    my_results_reduced = [file for file in os.listdir(output_dir) if file.endswith(".fasta") ]   
    for i in my_results_reduced:
        my_file = os.path.join(output_dir, i)
        if my_file != my_reference_IGRs:
            my_subjects.append(my_file)
    return my_subjects


def rbbh (evalue, threads, output_dir="results_sRNA/IGRs_blast"):
    print("Running Blast")
    my_reference_IGRs = my_genome(my_reference)
    my_subjects_IGRs = my_subjects(my_reference)
    
    #Forward BLASTn search
    for subject in my_subjects_IGRs:
        cline = NcbiblastnCommandline(query=my_reference_IGRs, 
                                  db=subject,
                                  evalue=evalue, 
                                  out = output_dir + "/" + my_reference_IGRs.split("/")[2] + "_vs_" + subject.split("/")[2]  + "_blastn_out.txt", 
                                  task = 'blastn',
                                  num_threads = threads,
                                  max_hsps = 1,
                                  max_target_seqs=1,
                                  #dust = 'no',
                                  outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue")
        stdout, stderr = cline()
    #Reverse  BLASTn search
    for subject in my_subjects_IGRs:
        cline = NcbiblastnCommandline(query=subject, 
                                  db=my_reference_IGRs,
                                  evalue=evalue, 
                                  out = output_dir + "/" + subject.split("/")[2] + "_vs_" + my_reference_IGRs.split("/")[2]  + "_blastn_out.txt", 
                                  task = 'blastn',
                                  num_threads = threads,
                                  max_hsps = 1,
                                  max_target_seqs=1,
                                  #dust = 'no',
                                  outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue")
        stdout, stderr = cline()
  
        
def analyze_blast_results(output_dir="results_sRNA/IGRs_blast"):
    rbbh_summary = []
    blast_result = [output for output in os.listdir(output_dir) if output.endswith(".txt") ] 
    my_reference_assembly = my_assembly_name(my_reference)

    forward_blast = []
    for i, j in enumerate(blast_result):
        if my_reference_assembly in blast_result[i].split("_vs_")[0]:
            forward_blast.append(j)

    reverse_blast = []
    for i, j in enumerate(blast_result):
        if my_reference_assembly not in blast_result[i].split("_vs_")[0]:
            reverse_blast.append(j)
    
    forward_blast = sorted(forward_blast)
    reverse_blast = sorted(reverse_blast)

    for forward, reverse in zip(forward_blast, reverse_blast):
        try:
            my_data_forward = pd.read_table(os.path.join(output_dir, forward), sep='\t', header=None)
            my_data_reverse = pd.read_table(os.path.join(output_dir, reverse), sep='\t', header=None)
        except EmptyDataError:
            print(f"Empty file found: {forward} or {reverse}. Skipping these files.")
            continue
        
        my_data_forward.columns = ["query", "subject", "identity", "coverage",
                                   "qlength", "slength", "alength",
                                   "bitscore", "E-value"]
        my_data_reverse.columns = ["query", "subject", "identity", "coverage",
                                   "qlength", "slength", "alength",
                                   "bitscore", "E-value"]
        
        my_data_forward['norm_bitscore'] = my_data_forward.bitscore / my_data_forward.qlength
        my_data_reverse['norm_bitscore'] = my_data_reverse.bitscore / my_data_reverse.qlength
        
        my_data_forward['qcov'] = my_data_forward.alength / my_data_forward.qlength
        my_data_reverse['qcov'] = my_data_reverse.alength / my_data_reverse.qlength
        my_data_forward['scov'] = my_data_forward.alength / my_data_forward.slength
        my_data_reverse['scov'] = my_data_reverse.alength / my_data_reverse.slength
        
        my_data_forward['qcov'] = my_data_forward['qcov'].clip(upper=1)
        my_data_reverse['qcov'] = my_data_reverse['qcov'].clip(upper=1)
        my_data_forward['scov'] = my_data_forward['scov'].clip(upper=1)
        my_data_reverse['scov'] = my_data_reverse['scov'].clip(upper=1)
        
        rbbh = pd.merge(my_data_forward, my_data_reverse[['query', 'subject']],
                        left_on='subject', right_on='query',
                        how='outer')
        
        rbbh = rbbh.loc[rbbh.query_x == rbbh.subject_y]
        rbbh = rbbh.groupby(['query_x', 'subject_x']).max()
        rbbh_summary.append(rbbh)
    
    rbbh_summary_coverage = []
    coverage = 0.7
    for rbbh in rbbh_summary:
        rbbh_coverage = rbbh[rbbh["qcov"] >= coverage]
        rbbh_summary_coverage.append(rbbh_coverage)
    
    rbbh_concat = pd.concat(rbbh_summary_coverage)
    return rbbh_concat

def concat_IGRs_reduced_fasta(output_dir = "results_sRNA/IGRs"):

    IGRs_fasta = [output for output in os.listdir(output_dir) if output.endswith(".fasta") ]     
    IGRs_database = os.path.join(output_dir, "IGRs_reduced_database.fasta")
    with open(IGRs_database , 'w') as outfile:
        for IGRs_sequences in IGRs_fasta:
            IGRs_file = os.path.join(output_dir, IGRs_sequences)
            with open(IGRs_file) as infile:
                for line in infile:
                    outfile.write(line)

#Extract the sequences conserved in the reference genome 
def extract_conserved_sequences_reference(my_reference, rbbh_dict, output_dir="results_sRNA/reference_results"):
    with open(os.path.join(output_dir, "headers_my_reference"), 'w') as output:
        for key, value in rbbh_dict.items():
            IGRs = value
            list_1 = IGRs.split(" ")
            list_2 = [key] + list_1
            if len(list_2) >= 4:
                output.write(str(list_2[0]) + '\n')
    
    input_file =  my_genome(my_reference)  # Aquí debes definir correctamente cómo obtener el archivo de entrada
    id_file = os.path.join(output_dir, "headers_my_reference")
    output_file = os.path.join(output_dir, "conserved_IGRs_my_reference.fasta")
    
    with open(id_file) as id_handle:
        wanted = set(line.rstrip("\n").split(None, 1)[0] for line in id_handle)
    print("Found %i unique identifiers in %s" % (len(wanted), id_file))
    
    records = (r for r in SeqIO.parse(input_file, "fasta") if r.description.split()[0] in wanted)
    count = SeqIO.write(records, output_file, "fasta")
    print("Saved %i records from %s to %s" % (count, input_file, output_file))
    if count < len(wanted):
        print("Warning %i IDs not found in %s" % (len(wanted) - count, input_file))

        
def clustalw_aln():
    headers = [header for header in os.listdir("results_sRNA/IGRs_cluster") if header.endswith(".headers") ]  

    for header in headers:
        input_file = os.path.join("results_sRNA/IGRs", "IGRs_reduced_database.fasta")
        id_file = os.path.join("results_sRNA/IGRs_cluster", header)
        output_file = id_file.split("headers")[0] + "fasta"
       
        with open(id_file) as id_handle:
            wanted = list(line.rstrip("\n").split(None, 1)[0] for line in id_handle)
        print("Found %i unique identifiers in %s" % (len(wanted), id_file))
    
        records = [r for r in SeqIO.parse(input_file, "fasta") if r.description in wanted]
        records_ordered= []
        for w in wanted:
            for j in records:               
                if j.id == w:
                    records_ordered.append(j)
        count = SeqIO.write(records_ordered, output_file, "fasta")
        print("Saved %i records from %s to %s" % (count, input_file, output_file))
        if count < len(wanted):
            print("Warning %i IDs not found in %s" % (len(wanted) - count, input_file))
    
        #muscle_cline = MuscleCommandline(muscle_exe, input=output_file, 
         #                                out=output_file.split()[0] + ".aln",
          #                               clwstrict=True)
        #stdout, stderr = muscle_cline()
        
        clustal_cline = ClustalwCommandline("clustalw", infile=output_file, 
                                         outfile=output_file.split()[0] + ".aln", 
                                         outorder="INPUT")
        stdout, stderr = clustal_cline()
        
        
        remove(id_file)
        

def my_assembly_name(my_reference):
    my_new_reference  = os.path.join(my_input_dir, my_reference)
    for seq_record in SeqIO.parse(my_new_reference, "genbank"):
        base_name=seq_record.dbxrefs[2]
        base_name_list=base_name.split(":")
        base_name_formated = base_name_list[1]
        assembly_name = os.path.splitext(os.path.basename(base_name_formated))[0]
    return assembly_name 

def coords(my_reference):  
    handle = my_input_dir+my_reference
    gene_locations = []    
    for seq_record in SeqIO.parse(handle, "genbank"):              
        for feature in seq_record.features:
            if feature.type == 'gene':
                my_start = feature.location.start.position
                my_end = feature.location.end.position
                gene_ID = str(feature.qualifiers['locus_tag'][0])
                gene_locations.append((gene_ID, my_start,my_end, seq_record.id))
                                                                
    outpath = "results_sRNA/TransTermHP/"+  my_assembly_name(my_reference) + ".coords"
    with open(outpath, 'w') as output:
        for row in gene_locations:
                    output.write(str(row[0])+chr(9)+str(row[1])+chr(9)+str(row[2])+chr(9)+str(row[3])+'\n')


def fake_coords(handle, strand):
    fake = []
    for seq_record in SeqIO.parse(handle, "fasta"):
        
        data = {"fakegene1":[1, 2, seq_record.id], "fakegene2":[len(seq_record)-1, len(seq_record),  seq_record.id]}
        fake_genes = pd.DataFrame.from_dict(data, orient='index')
        fake.append(fake_genes)
    
    fake_concat = pd.concat(fake)
    fake_concat.to_csv ("results_sRNA/TransTermHP/"+ my_assembly_name(my_reference) + "_"+ strand + "_fake.coords", index = True, header = False, sep="\t")


def transtermhp(fasta, coords):
    cmd = ["transterm", '-p', 'expterm.dat', fasta, coords]
    output = subprocess.check_output(cmd, text=True)
    #   TERM 1         4342 - 4366     + F    93 -11.5 -3.22878 | opp_overlap 4342, overlap 4340 4357
    ttre = re.compile(
        r'^  (?P<name>.*) (?P<start>\d+) - (?P<end>\d+)\s+'
        r'(?P<strand>[-+])\s+(?P<loc>[GFRTHNgfr]+)\s+'
        r'(?P<conf>\d+)\s+(?P<hp>[0-9.-]+)\s+(?P<tail>[0-9.-]+)'
    )
    
    features = [] 
    record_fasta_replicon = SeqIO.read(fasta, "fasta")        
    batches = output.split('SEQUENCE ')
    batch_id = record_fasta_replicon.id
    for batch in batches[1:]:
        batch_lines = batch.split('\n')
        # Strip the header
        interesting = batch_lines[2:]
        unformatted = [x for x in interesting if x.startswith('  ')][0::2]
        for terminator in unformatted:
            m = ttre.match(terminator)
            if m:
                start = int(m.group('start')) - 1
                end = int(m.group('end'))
                if m.group('strand') == '+':
                    strand = 1
                else:
                    strand = 0
                    
                feature = {"Sequence":batch_id, "ID": m.group('name'), 
                           "Strand":strand, "Start":start, "End":end, 
                                "type":"terminator", "score": m.group('conf'),                           
                                "source": "TransTermHP_2.09"} 
                features.append(feature)
                
    df = pd.DataFrame(features)
    return df

def promoters_in_IGRs(fasta_file, promoters, output):
    headers_IGRs_conserved = []
    flanking_genes = []
    coordinates = []
    orientation  = []
    
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        headers_IGRs_conserved.append(seq_record.id)
        flanking_genes.append(seq_record.id.split("++")[1])
        coordinates.append(seq_record.id.split("++")[2])
        orientation.append(seq_record.id.split("++")[3])
    
    IGRs_conserved_dict2 = dict(zip(headers_IGRs_conserved, coordinates))
    
    # Locate promoters predicted by G4PromFinder or Promotech in the conserved IGRs
    empty_lists = [[] for _ in range(len(coordinates))]
    promoters_in_IGRs = {headers_IGRs_conserved[i]: empty_lists[i] for i in range(len(headers_IGRs_conserved))}
    
    print("Locating promoters in IGRs sequences")
    
    for key, value in IGRs_conserved_dict2.items():
        start_IGR = int(IGRs_conserved_dict2[key].split("-")[0])
        end_IGR = int(IGRs_conserved_dict2[key].split("-")[1])
        
        for j in promoters.itertuples():
            j_start = int(j.start)
            j_end = int(j.end)
            
            if start_IGR <= j_start and end_IGR >= j_end and key.split("++")[0] == j.chrom:
                promoter_region = (j.start, j.end, j.strand, j.sequence)
                promoters_in_IGRs[key].append(promoter_region)
        
    """without_promoter = []
    for key, value in promoters_in_IGRs.items():
        if len(value) == 0:
            without_promoter.append(1)
            
    number_of_promoters = []
    for key, value in promoters_in_IGRs.items():
        number_of_promoters.append(len(value))"""
    
    # Save promoters in IGRs
    df2 = pd.DataFrame.from_dict(promoters_in_IGRs, orient="index")
    new_names = {str(col): f'Promoter{int(col)+1}' for col in df2.columns}
    new_names_list = list(new_names.values())
    df2.columns = new_names_list
    df2.to_csv(os.path.join("results_sRNA/reference_results", output), index=True, header=False)
    
    return df2


def terminator_in_IGRs(fasta_file, terminators, output):
    headers_IGRs_conserved = []
    flanking_genes = []
    coordinates = []
    orientation  = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        headers_IGRs_conserved.append(seq_record.id)
        flanking_genes.append(seq_record.id.split("++")[1])
        coordinates.append(seq_record.id.split("++")[2])
        orientation.append(seq_record.id.split("++")[3])
    
    IGRs_conserved_dict2 = dict(zip(headers_IGRs_conserved, coordinates))




    empty_lists2 =  [[] for _ in range(len(coordinates))]
    terminators_in_IGRs = {headers_IGRs_conserved[i]:empty_lists2[i] for i in range(len(headers_IGRs_conserved))}

    print("Locating terminators in IGRs sequences")

    for key, value in IGRs_conserved_dict2.items():
        start_IGR = int(IGRs_conserved_dict2[key].split("-")[0])
        end_IGR = int(IGRs_conserved_dict2[key].split("-")[1])
        for j in terminators_transtermhp.itertuples():
            if start_IGR <= j.Start and end_IGR >= j.End and key.split("++")[0] == j.Sequence:
                terminator_region = (j.Start, j.End, j.Strand)
                terminators_in_IGRs[key].append(terminator_region)

    #Save terminators in IGRs
    df3= pd.DataFrame.from_dict(terminators_in_IGRs, orient="index")
    new_names = {str(col): f'Terminator{int(col)+1}' for col in df3.columns}
    new_names_list = list(new_names.values())
    df3.columns = new_names_list


    df3.to_csv(output, index = True, header = False)

    #Check the number of conserved IGRs with terminator 
    """has_terminators= []
    for key, value in terminators_in_IGRs.items():
        if len(value) != 0:
            has_terminators.append(value)

    without_terminator = []
    for key, value in terminators_in_IGRs.items():
        if len(value)==0:
            without_terminator.append(1)

    with_terminator = len(terminators_in_IGRs) - len(without_terminator)"""

    return df3


def shuffle_sequences(fasta_file):
    shuffled_sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Shuffle the sequence
        shuffled_seq = "".join(random.sample(str(record.seq), len(record.seq)))
        # Create a new SeqRecord with the shuffled sequence
        shuffled_record = record[:]
        shuffled_record.id =  record.id + "_shuffled"
        shuffled_record.description =  record.description + "_shuffled"
        shuffled_record.seq = shuffled_seq
        shuffled_sequences.append(shuffled_record)
    return shuffled_sequences
                    



    
#%%  Workflow


# Define a function to log directory contents
def log_directory_contents(directory):
    print(f"Contents of directory '{directory}':")
    for filename in os.listdir(directory):
        print(filename)

def workflow(my_input_dir, my_reference):
  
    # Define directories and create them if necessary
    my_current_dirs = [name for name in os.listdir(".") if os.path.isdir(name)]
    
    if "results_sRNA" in my_current_dirs:
        shutil.rmtree("results_sRNA")
    os.mkdir("results_sRNA")
    
    
    
    os.mkdir(os.path.join("results_sRNA", "IGRs_blast"))
    os.mkdir(os.path.join("results_sRNA", "IGRs_cluster"))
    os.mkdir(os.path.join("results_sRNA", "IGRs_alignments"))
    os.mkdir(os.path.join("results_sRNA", "reference_results"))
    os.mkdir(os.path.join("results_sRNA", "IGRs"))
    
    """   
     # Log the contents of the input directory
        log_directory_contents(my_input_dir)
    
        # Check if the reference file exists
        reference_file_path = os.path.join(my_input_dir, my_reference)
        if not os.path.exists(reference_file_path):
            print(f"Reference file '{my_reference}' not found in directory '{my_input_dir}'")
            return
    """
        # Extract IGRs for the genomes in the input directory
    my_files = os.listdir(my_input_dir)
    for i in my_files:
        print(f"Processing file {i}")
        get_intergenic_regions(os.path.join(my_input_dir, i), intergenic_length=50, max_intergenic_length=600)
    
        # Summarize the number of IGRs in each genome
    IGRs_per_genome = []
    my_results = os.listdir(os.path.join("results_sRNA", "IGRs"))
    for i in my_results:
        num = len([1 for line in open(os.path.join("results_sRNA", "IGRs", i)) if line.startswith(">")])
        print(num, "IGRs in file", i)
        IGRs_per_genome.append(num)
    
        # Cut the borders of the sequence to avoid including regulatory elements
    for i in my_results:
        print(f"Processing file {i}")
        trim_sequences(os.path.join("results_sRNA", "IGRs", i), 0, 0)
        remove(os.path.join("results_sRNA", "IGRs", i))



    fasta_file  = my_genome(my_reference, output_dir = "results_sRNA/IGRs")
    shuffled_sequences = shuffle_sequences(fasta_file)
    output_file = fasta_file
    SeqIO.write(shuffled_sequences, output_file, "fasta")



    # Create the databases
    make_blast_db()

    # Run reciprocal best BLAST
    threads = 18
    evalue = 1e-5
    rbbh(evalue, threads)

    # Analyze BLAST results
    rbbh_concat = analyze_blast_results()

    # Save RBBH results
    rbbh_concat.to_csv('results_sRNA/IGRs_blast/RBBH.csv', index=True, header=True)

    # Get the RBBH ids for multiple alignment in subsequent steps
    rbbh_concat.reset_index(inplace=True)
    rbbh_id = rbbh_concat[['query_x', 'subject_x']]
    rbbh_id = rbbh_id.groupby('query_x')['subject_x'].agg(' '.join)
    rbbh_dict = rbbh_id.to_dict()

    # Extract the sequences conserved in the reference genome
    extract_conserved_sequences_reference(my_reference, rbbh_dict)
    
    #Create the headers for the sequences conserved in IGR clusters
    for key, value in rbbh_dict.items():
        IGRs =  value
        list_1 = IGRs.split(" ")
        list_2 = [key] + list_1
        file =  os.path.join("results_sRNA/IGRs_cluster", key + ".headers")
        
        if len(list_2) >=4:
            with open(file, 'w') as output:
                for row in list_2:
                    output.write(str(row) + '\n')
                    
    #Extract IGR clusters and align them

    concat_IGRs_reduced_fasta() 
    clustalw_aln()
    
    #Move the alignments to the "IGRs_alignments" directory
    alignments = [aln for aln in os.listdir("results_sRNA/IGRs_cluster") if aln.endswith(".aln")]  
    for file_name in alignments:
        full_file_name = os.path.join("results_sRNA/IGRs_cluster", file_name)
        shutil.move(full_file_name, os.getcwd())




#%%

# Define la función workflow y otras funciones necesarias

# Define los parámetros
my_input_dir = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/Genbank_clade2_pangenome"
my_reference = "GCF_005519465.1_ASM551946v1_genomic.gbff"


# Define el número de veces que quieres ejecutar la función
num_executions = 1000000  # Cambia este valor según lo que necesites
start_time = time.time()
for _ in range(num_executions):
    workflow(my_input_dir, my_reference)

elapsed_time = time.time() - start_time

# Imprime el tiempo transcurrido
print(f"Tiempo total de ejecución: {elapsed_time} segundos")
