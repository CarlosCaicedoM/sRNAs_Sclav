#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 22:45:32 2022

@author: Carlos Caicedo-Montoya
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
import wget
import gzip
from Bio.SeqUtils import gc_fraction
from collections import Counter
from scipy.stats import chi2_contingency
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
    for i,j in enumerate(blast_result):
        if my_reference_assembly in blast_result[i].split("_vs_")[0]:
            forward_blast.append(j)

    reverse_blast = []
    for i,j in enumerate(blast_result):
        if my_reference_assembly not in blast_result[i].split("_vs_")[0]:
            reverse_blast.append(j)
    
    forward_blast = sorted(forward_blast)
    reverse_blast = sorted(reverse_blast)

    for forward, reverse in zip(forward_blast, reverse_blast):
        # The detailed explanation for this chunk of code can be found at
        # https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/05-blast_for_rbh.html
        
        my_data_forward = pd.read_table(os.path.join(output_dir, forward), sep='\t', header=None)
        
        my_data_reverse = pd.read_table(os.path.join(output_dir, reverse), sep='\t', header=None)
        
        
        my_data_forward.columns = ["query", "subject", "identity", "coverage",
                   "qlength", "slength", "alength",
                   "bitscore", "E-value"]
        
        
        my_data_reverse.columns = ["query", "subject", "identity", "coverage",
                   "qlength", "slength", "alength",
                   "bitscore", "E-value"]
        
        my_data_forward['norm_bitscore'] = my_data_forward.bitscore/my_data_forward.qlength
        my_data_reverse['norm_bitscore'] = my_data_reverse.bitscore/my_data_reverse.qlength
        
        
        # Create query and subject coverage columns in both dataframes
        my_data_forward['qcov'] = my_data_forward.alength/my_data_forward.qlength
        my_data_reverse['qcov'] = my_data_reverse.alength/my_data_reverse.qlength
        my_data_forward['scov'] = my_data_forward.alength/my_data_forward.slength
        my_data_reverse['scov'] = my_data_reverse.alength/my_data_reverse.slength
        
        # Clip maximum coverage values at 1.0
        my_data_forward['qcov'] = my_data_forward['qcov'].clip(upper = 1)
        my_data_reverse['qcov'] = my_data_reverse['qcov'].clip(upper = 1)
        my_data_forward['scov'] = my_data_forward['scov'].clip(upper = 1)
        my_data_reverse['scov'] = my_data_reverse['scov'].clip(upper = 1)
        
        
        # Merge forward and reverse results
        rbbh = pd.merge(my_data_forward, my_data_reverse[['query', 'subject']],
                        left_on='subject', right_on='query',
                        how='outer')
        
        rbbh = rbbh.loc[rbbh.query_x == rbbh.subject_y]
        
        # Group duplicate RBH rows, taking the maximum value in each column
        rbbh = rbbh.groupby(['query_x', 'subject_x']).max()
    
        #Append the results to RBBH_summary
        rbbh_summary.append(rbbh)
    
    #Filter by Coverage 
    rbbh_summary_coverage =[]
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
def extract_conserved_sequences_reference(output_dir = "results_sRNA/reference_results"):
    with open(os.path.join(output_dir, "headers_my_reference"), 'w') as output:
        for key, value in rbbh_dict.items():
            IGRs =  value
            list_1 = IGRs.split(" ")
            list_2 = [key] + list_1
            if len(list_2) >=4:
                output.write(str(list_2[0]) + '\n')
            
    input_file = my_genome(my_reference)
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
                    
#%% Find promoters with G4PromFinder

def G4PromFinder(File_genome):    
    try:
        my_reference_G4PromFinder = "results_sRNA/G4PromFinder/" + my_assembly_name(my_reference) + "_Promoter-coordinates.txt"
        
        q = open(my_reference_G4PromFinder ,"a")
        q.write("Sequence,Region,Strand,Start,End")
        for seq_record in SeqIO.parse(File_genome, "fasta"):
            print(len(seq_record))
            mc = []    
            mc.append(str(seq_record.seq))
            print ("G4PromFinder is working on sequence %s, please wait..." % (seq_record.id))
            genome = mc[0] 
            genome = genome.upper()
            gc_total = 100*(genome.count("G")+genome.count("C"))/len(genome)
            pattern1 = "G{2,4}\D{1,10}G{2,4}\D{1,10}G{2,4}\D{1,10}G{2,4}"
            pattern2 = "C{2,4}\D{1,10}C{2,4}\D{1,10}C{2,4}\D{1,10}C{2,4}"
            patternA = "TA\D{3}T"
            patternB = "A\D{3}TA"
            patternC = "TTGAC"
            patternD = "GTCAA"
            predictions = 0
            cod1 = 1
            cod2 = 1
            AT_content = []
            GC_content = []
            SI1 = 0
            SI2 = 0
            NO1 = 0
            NO2 = 0
            #q = open("Promoter_coordinates1.txt","a")
            #q.write("\nRegion,Strand,Start,End")
            
            scale = 50
            for j in range(100000000):
                x1 = scale + j
                x2 = scale + j + 25
                OK = 0
                if x2 > len(genome):
                    break
                else:
                    w = genome[x1:x2]
                    at =100*(w.count("A")+w.count("T"))/len(w)
                    prom = genome[x1-50:x2]
                    if at >= 40:
                        if re.findall(pattern1,prom) or re.findall(pattern2,prom):
                            if re.findall(pattern1,prom):
                                ricerca = re.findall(pattern1,prom)
                                for lk in ricerca:
                                    if len(lk) <= 30:
                                        OK += 1
                                        break
                            elif re.findall(pattern2,prom):
                                ricerca = re.findall(pattern2,prom)
                                for lk in ricerca:
                                    if len(lk) <= 30:
                                        OK += 1
                                        break
                            if OK == 1:
                                predictions += 1
                                dati = []
                                for g in range(50):
                                    if x2+g <= len(genome):
                                        w = genome[x1+g:x2+g]
                                        at =100*(w.count("A")+w.count("T"))/len(w)
                                        prom = genome[x1+g-50:x2+g]
                                        if at >= 40:
                                            if re.search(pattern2,prom) or re.search(pattern1,prom):
                                                dati.append(at)
                                            else:
                                                dati.append(0)
                                        else:
                                            dati.append(0)
                                    else:
                                        break
                                maxP = np.argmax(dati)
                                at = np.max(dati)
                                x1 = x1 + maxP
                                x2 = x2 + maxP
                                scale = x2 - j - 1 + 50
                                prom = genome[x1-50:x2]
                                gc = 100*(prom.count("G")+prom.count("C"))/len(prom)
                                q.write("\n")
                                q.write(seq_record.id)
                                q.write(",")
                                q.write("P")
                                q.write(str(cod1))
                                q.write("plus,")
                                q.write("positivo,")
                                q.write(str(x1-50))
                                q.write(",")
                                q.write(str(x2))
                                cod1 += 1
                                AT_content.append(at)
                                GC_content.append(gc)
                                if re.search(patternA,prom):
                                    SI1 += 1
                                else:
                                    NO1 += 1
                                if re.search(patternC,prom):
                                    SI2 += 1
                                else:
                                    NO2 += 1
                                    
            scale = 0
            for j in range(100000000):
                x1 = scale + j
                x2 = scale + j + 25
                OK = 0
                if x2 > len(genome)-50:
                    break
                else:
                    w = genome[x1:x2]
                    at =100*(w.count("A")+w.count("T"))/len(w)
                    prom = genome[x1:x2+50]
                    if at >= 40:
                        if re.findall(pattern1,prom) or re.findall(pattern2,prom):
                            if re.findall(pattern1,prom):
                                ricerca = re.findall(pattern1,prom)
                                for lk in ricerca:
                                    if len(lk) <= 30:
                                        OK += 1
                                        break
                            elif re.findall(pattern2,prom):
                                ricerca = re.findall(pattern2,prom)
                                for lk in ricerca:
                                    if len(lk) <= 30:
                                        OK += 1
                                        break
                            if OK == 1:
                                dati = []
                                predictions += 1
                                for g in range(50):
                                    if x2+g <= len(genome):
                                        w = genome[x1+g:x2+g]
                                        at =100*(w.count("A")+w.count("T"))/len(w)
                                        prom = genome[x1+g:x2+g+50]
                                        if at >= 40:
                                            if re.search(pattern1,prom) or re.search(pattern2,prom):
                                                dati.append(at)
                                            else:
                                                dati.append(0)
                                        else:
                                            dati.append(0)
                                    else:
                                        break
                                maxP = np.argmax(dati)
                                at = np.max(dati)
                                x1 = x1 + maxP
                                x2 = x2 + maxP 
                                scale = x2 - j - 1 + 50
                                prom = genome[x1:x2+50]
                                gc = 100*(prom.count("G")+prom.count("C"))/len(prom)
                                q.write("\n")
                                q.write(seq_record.id)
                                q.write(",")
                                q.write("P")
                                q.write(str(cod2))
                                q.write("minus,")
                                q.write("negativo,")
                                q.write(str(x1))
                                q.write(",")
                                q.write(str(x2+50))
                                cod2 += 1
                                AT_content.append(at)
                                GC_content.append(gc)
                                if re.search(patternB,prom):
                                    SI1 += 1
                                else:
                                    NO1 += 1
                                if re.search(patternD,prom):
                                    SI2 += 1
                                else:
                                    NO2 += 1
                                    
           
            
            with open("results_sRNA/G4PromFinder/" + my_assembly_name(my_reference) + "_ABOUT_PROMOTERS.txt", "a") as file:
                file.write("Sequence:")
                file.write(seq_record.id)
                file.write("\nMean GC% content of sequence:")
                file.write(str(gc_total))
                file.write("%")
                file.write("\nTotal number of prediction in the analyzed sequence:")
                file.write(str(predictions))
                file.write("\nMean GC-content of putative promoters:")
                file.write(str(np.mean(GC_content)))
                file.write("%\nMean AT-content of the putative promoters AT-rich elements:")
                file.write(str(np.mean(AT_content)))
                file.write("%\nPercent of putative promoters with the -10 consensus TANNNT:")
                file.write(str(100*SI1/(SI1+NO1)))
                file.write("%\nPercent of putative promoters with the -35 consensus TTGAC:")
                file.write(str(100*SI2/(SI2+NO2)))
                file.write("%\n")
        q.close()
        print("Work finished, see output files in the directory results_sRNA/G4PromFinder/")
        
    except IOError:   
        print ("File %s inexistent in the current directory!" %(File_genome))


#%% Gbk_to_ppt

"""
Title: GenBank file to Protein Table file (.ptt) and RNA table file (.rnt) parser
Input files: [./sequence.gb, ./*.gb (batch)]

This script creates a Protein Table file (.ptt) and a RNA table file (.rnt) from the given GenBank file
Multiple files can be given (or using *) for batch processing

@author: Henry Wiersma UMCG, Groningen
@date: 7-11-2017
@version: 1.0
"""

#createFile(pttOutputFilePath, pttHeader, description, length, pttRows)
#filePath = pttOutputFilePath
#headerTemplate=pttHeader
#rows = pttRows

def createFile(filePath, headerTemplate, description, length, rows):
    outputFile = open(filePath, 'w')
    header = headerTemplate.format(description=description,
                                   length=length,
                                   numRows=len(rows))
    outputFile.write(header)
    for row in rows:
        outputFile.write("%s\n" % row)
    outputFile.close()


def processFile(filePath, outputPath):
    #rna feature types in the GenBank file (https://www.ncbi.nlm.nih.gov/books/NBK293913/)
    rntFeatureTypes = ["rna", "mRNA", "tRNA", "rRNA", "ncRNA", "tmRNA", "misc_RNA"]

    #Template of the header of the tables
    pttHeader="{description} - 1..{length}\n{numRows} proteins\nLocation\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n"
    rntHeader="{description} - 1..{length}\n{numRows} RNAs\nLocation\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n"

    #template of the rows of the table
    #row="{location}\t{strand}\t{length}\t{pid}\t{gene}\t{synonym}\t{code}\t{cog}\t{product}\r\n"
    row="{location}\t{strand}\t{length}\t{pid}\t{gene}\t{synonym}\t{code}\t{cog}\t{product}"


    file = SeqIO.parse(filePath, "genbank")
    
    print("Process file", os.path.basename(filePath))
    for seqRec in file:
        fileNameBase = seqRec.id
        pttRows = []
        rntRows = []

        # get description of the gb file
        description = "Description not available"
        if seqRec.description and seqRec.description != "" and seqRec.id and seqRec.id != "":
            description = seqRec.id + " " + seqRec.description
        elif seqRec.description and seqRec.description != "":
            description = seqRec.description
        elif seqRec.id and seqRec.id != "":
            description = seqRec.id

        # get the sequence length
        length = 0
        if seqRec.seq:
            length = len(seqRec.seq)
        a =  randint(100000000, 500000000)     # randint is inclusive at both ends
        pid = a
        for i, feature in enumerate(seqRec.features):
            # location
            paramLocation = "{}..{}".format((feature.location.start + 1), feature.location.end)

            # strand
            paramStrand = "+" if feature.location.strand == 1 else "-"

            # gene name
            paramGene = "-"
            if "gene" in feature.qualifiers and len(feature.qualifiers["gene"]) > 0 and feature.qualifiers["gene"][
                0] != "":
                paramGene = (feature.qualifiers["gene"][0])

            # locus tag (Synonym)
            paramSynonym = "-"
            if "locus_tag" in feature.qualifiers and len(feature.qualifiers["locus_tag"]) > 0 and \
                            feature.qualifiers["locus_tag"][0] != "":
                paramSynonym = (feature.qualifiers["locus_tag"][0])

            # product descriptopn
            paramProduct = "-"
            if "product" in feature.qualifiers and len(feature.qualifiers["product"]) > 0 and \
                            feature.qualifiers["product"][
                                0] != "":
                paramProduct = (feature.qualifiers["product"][0])
                
            if feature.type == "CDS":

                # length of product
                paramLength = round((feature.location.end - feature.location.start) / 3 - 1)

                tempRow = row.format(location=paramLocation,
                                     strand=paramStrand,
                                     length=paramLength,
                                     pid=pid+ int(i/2),
                                     gene=paramGene,
                                     synonym=paramSynonym,
                                     code="-",
                                     cog="-",
                                     product=paramProduct)
                pttRows.append(tempRow)

            elif feature.type in rntFeatureTypes:

                # length of product
                paramLength = round((feature.location.end - feature.location.start) - 1)

                tempRow = row.format(location=paramLocation,
                                     strand=paramStrand,
                                     length=paramLength,
                                     pid="-",
                                     gene=paramGene,
                                     synonym=paramSynonym,
                                     code="-",
                                     cog="-",
                                     product=paramProduct)
                rntRows.append(tempRow)

        # create files
       
        pttOutputFilePath = os.path.join(outputPath, "{}.ptt".format(fileNameBase))
        rntOutputFilePath = os.path.join(outputPath, "{}.rnt".format(fileNameBase))

        print("number of ptt rows:\t{}\t({})".format(len(pttRows), pttOutputFilePath))
        print("number of rnt rows:\t{}\t({})".format(len(rntRows), rntOutputFilePath))

        createFile(pttOutputFilePath, pttHeader, description, length, pttRows)
        createFile(rntOutputFilePath, rntHeader, description, length, rntRows)




    
#%%  Workflow

start_time = time.time()

##Parameters for the main function
#Define the folder with the genomes to analyze

#my_input_dir = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/Genbank_group1"
my_input_dir = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/Genbank_clade2_pangenome"
#my_input_dir = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/Genbank_clade2_pangenome_reduced"



my_reference = "GCF_005519465.1_ASM551946v1_genomic.gbff"

#Extract IGRs from the genomes to compare
#Define a folder to store temporary files and intermediate results

my_current_dirs = [name for name in os.listdir(".") if os.path.isdir(name)]

if "results_sRNA" in my_current_dirs:
    shutil.rmtree("results_sRNA")
    os.mkdir("results_sRNA")
else:
    os.mkdir("results_sRNA")

os.mkdir(os.path.join("results_sRNA", "IGRs_blast"))
os.mkdir(os.path.join("results_sRNA", "IGRs_cluster"))
os.mkdir(os.path.join("results_sRNA", "IGRs_alignments"))
os.mkdir(os.path.join("results_sRNA", "reference_results"))  
os.mkdir(os.path.join("results_sRNA", "IGRs"))

#to get the data consider to use the following comand
"""
pip install ncbi-genome-download
ncbi-genome-download bacteria -F genbank -A Streptomyces_Assemby_Accession_Group1_Sclav.txt -v

#To unzip multiple files
chmod +x unzip_several.sh
./unzip_several.sh
"""

#Extract IGRs for the genomes in my input directory
my_files = os.listdir(my_input_dir)
for i in my_files:
    print("Processing file %s" %i )
    get_intergenic_regions(os.path.join(my_input_dir, i), intergenic_length=50, max_intergenic_length=600)

#Summarize the number of IGRs in each genome
IGRs_per_genome = []
my_results = os.listdir(os.path.join("results_sRNA", "IGRs"))
for i in my_results:
    num = len([1 for line in open(os.path.join("results_sRNA", "IGRs", i)) if line.startswith(">")])    
    print(num, "IGRs in file", i)                                
    IGRs_per_genome.append(num)

#Cut the borders of the sequence to avoid to include regulatory elements in the analysis
#that can belong to the 5' or 3' UTRs
for i in my_results:
   print("Processing file %s" %i )
   trim_sequences(os.path.join("results_sRNA", "IGRs", i), 0, 0)
   remove(os.path.join("results_sRNA", "IGRs", i))


## Create the databases
make_blast_db()

         
### Run reciprocal best blast 
threads = 18
evalue = 1e-5
rbbh(evalue, threads)


### Analyze_blast_results
rbbh_concat=analyze_blast_results()   
    
#Save RBBH results
rbbh_concat.to_csv ('results_sRNA/IGRs_blast/RBBH.csv', index = True, header=True)
    
    
#Get the RBBH ids in order to carry out multiple alignment in subsequent steps  
rbbh_concat.reset_index(inplace=True)  
rbbh_id = rbbh_concat[['query_x', 'subject_x']]  
rbbh_id = rbbh_id.groupby('query_x')['subject_x'].agg(' '.join)
rbbh_dict = rbbh_id.to_dict()
        
#Extract the sequences conserved in the reference genome 

extract_conserved_sequences_reference()         


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
alignments = [aln for aln in os.listdir("results_sRNA/IGRs_cluster") if aln.endswith(".aln") ]  
for file_name in alignments:
    full_file_name = os.path.join("results_sRNA/IGRs_cluster", file_name)
    shutil.move(full_file_name, "results_sRNA/IGRs_alignments")


# Convert genbank to fasta for the reference genome
my_reference_fasta = "results_sRNA/reference_results/" + my_assembly_name(my_reference) +".fasta"
records = SeqIO.parse(os.path.join(my_input_dir,  my_reference),  "genbank")
count = SeqIO.write(records, my_reference_fasta , "fasta")
print("Converted %i records" % count)


#  Find promoters with G4PromFinder
os.mkdir(os.path.join("results_sRNA","G4PromFinder"))
G4PromFinder(my_reference_fasta)


#Call promoter results from G4PromFinder or promotech
promoters = pd.read_table("results_sRNA/G4PromFinder/" + my_assembly_name(my_reference)+ "_Promoter-coordinates.txt", sep=",")
 
promotech_results = "Promotech"
promotech = [prom for prom in os.listdir(promotech_results) if prom.endswith(".csv") ]  

promotech_table = []
for i in promotech:
    predictions = pd.read_csv("Promotech/" + i, sep="\t")
    promotech_table.append(predictions)

promoters_promotech = pd.concat(promotech_table, axis=0)


#Read the data from the conserved IGRs
IGRs_conserved_reference = os.path.join("results_sRNA", "reference_results", "conserved_IGRs_my_reference.fasta")
#IGRs_conserved_reference = os.path.join("results_sRNA_clade2_pangenome_gene_nontrimmed_50-600", 
#                                       "reference_results", "conserved_IGRs_my_reference.fasta")

promoters_in_IGRs_conserved= promoters_in_IGRs(IGRs_conserved_reference, 
                                      promoters_promotech,
                                      "Promoters_in_conserved_IGRs")


#Determine terminators in the genome
os.mkdir(os.path.join("results_sRNA", "TransTermHP"))

#This function was inspired by the code found here (although modified for the purposes)
#of the current approach
#https://github.com/galaxyproject/tools-iuc/blob/master/tools/transtermhp/transtermhp.py


#convert genbank format to ptt
file = os.path.join(my_input_dir, my_reference)
outputPath = os.path.join("results_sRNA", "TransTermHP")
processFile(file, outputPath)

#Extract fasta sequence of each replicon
for seq_req in SeqIO.parse(file, "genbank"):
    print(seq_req.id)
    SeqIO.write(seq_req, "results_sRNA/reference_results/"+ seq_req.id + ".fasta", "fasta")


#run transtermhp on each replicon
replicons = []
for seq_req in SeqIO.parse(file, "genbank"):
    replicons.append(seq_req.id)

transterm_table = []
for replicon in replicons:
    coords_ptt = os.path.join(outputPath, replicon) + ".ptt"
    fasta = os.path.join("results_sRNA/reference_results/", replicon) + ".fasta"
    transterm_ptt=transtermhp(fasta, coords_ptt)
    transterm_table.append(transterm_ptt)
terminators_transtermhp = pd.concat(transterm_table, axis=0)

#Locate terminators predicted by TrasnTermHP in the conserved IGRs
terminator_results_conserved_IGRS_path = os.join.path("results_sRNA", 
                                                      "TransTermHP", 
                                                      "terminator_in_conserved_IGRs")

terminator_in_conserved_IGRs = terminator_in_IGRs(IGRs_conserved_reference, 
                                        terminators_transtermhp,                                         
                                        terminator_results_conserved_IGRS_path)



#Run RNAz using the conserved IGR sequences
#subprocess.call('./run_RNAz_local.sh')
subprocess.call('./run_RNAz_local_parallel.sh')


rnaz_results = [rnaz for rnaz in os.listdir("results_sRNA/RNAz_out") if rnaz.endswith(".out") ]

RNAz_all = os.path.join("results_sRNA/RNAz_out", "RNAz_all.out")
with open(RNAz_all , 'w') as outfile:
    for rnaz in rnaz_results:
        RNAz_file = os.path.join("results_sRNA/RNAz_out", rnaz)
        with open(RNAz_file) as infile:
            for line in infile:
                outfile.write(line)

#Extract identifiers of reference IGRs with conserved secondary structure 
sequences_reference_rnaz= set(line for line in open(RNAz_all) if (line.startswith(">NZ_CP027858.1") or line.startswith(">NZ_CP027859.1")))
sequences_reference_rnaz_list = list(sequences_reference_rnaz)
temp=[l.strip('\n\r') for l in sequences_reference_rnaz_list ]


#Get FASTA sequence of IGRs with conserved secondary structure 
input_file = os.path.join("results_sRNA/reference_results", "conserved_IGRs_my_reference.fasta")
output_file = os.path.join("results_sRNA/reference_results", "conserved_IGRs_my_reference_RNAz.fasta")    
records_rnaz = [r for r in SeqIO.parse(input_file, "fasta") if (">"+r.description.split("-")[0]+"-C") in temp]
count = SeqIO.write(records_rnaz, output_file, "fasta")
print("Saved %i records from %s to %s" % (count, input_file, output_file))
if count < len(temp):
    print("Warning %i IDs not found in %s" % (len(input_file) - count, input_file))


#Search IGRs with conserved secondary structure in the RFAM database
my_conserved_IGRs_RNAz = os.path.join("results_sRNA/reference_results", "conserved_IGRs_my_reference_RNAz.fasta")
rfam_db = os.path.join(os.getcwd(), "raw_data", "rfam_sequences", "Rfam.fa")

#BLASTn search
cline_rfam = NcbiblastnCommandline(query=my_conserved_IGRs_RNAz, 
                                  db=rfam_db,
                                  evalue=1e-6, 
                                  out = os.path.join("results_sRNA", "reference_results", "my_conserved_IGRs_RNAz_blastn_out.txt"), 
                                  task = 'blastn',
                                  num_threads = 10,
                                  max_hsps = 1,
                                  max_target_seqs=1,
                                  #dust = 'no',
                                  outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue")
stdout, stderr = cline_rfam()


rfam_results = pd.read_table(os.path.join("results_sRNA", "reference_results", "my_conserved_IGRs_RNAz_blastn_out.txt"), sep='\t', header=None)

rfam_results.columns = ["query", "subject", "identity", "coverage",
                   "qlength", "slength", "alength",
                   "bitscore", "E-value"]
rfam_results.set_index("query", inplace=True)

#concatenate the results of Promotech and TranstermHP in the conserved IGRs

prom_and_term = pd.concat([promoters_in_IGRs_conserved, 
                           terminator_in_conserved_IGRs[0]], axis = 1)


rnaz_ids = []
for rnaz in records_rnaz:
    rnaz_ids.append(rnaz.id)

have_rnaz =  [True for _ in range(len(rnaz_ids))]

rnaz_dict = dict(zip(rnaz_ids, have_rnaz))
rnaz_df = pd.DataFrame.from_dict(rnaz_dict, orient="index")
rnaz_df.columns = ["RNAz"]


prom_term_RNAz = pd.concat([prom_and_term, rnaz_df], axis = 1)


putative_sRNAs_rfam_blast = pd.concat([prom_term_RNAz, rfam_results], axis = 1)


#Annotate with RFAM
#url = 'https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz'
#filename = wget.download(url)
#url2="https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin"
#filename2  = wget.download(url2)

#Decompress covariance models
#with gzip.open(filename, 'rb') as f_in:
 #   with open(' Rfam.cm', 'wb') as f_out:
  #      shutil.copyfileobj(f_in, f_out)


os.mkdir("results_sRNA/Infernal")
Infernal_Results = os.path.join("results_sRNA/Infernal","Infernal_Sclav")
#cmd_cmpress = ["cmpress", 'Rfam_2024-04.cm']
#subprocess.call(cmd_cmpress)

# Determine the total database size for the genome you are annotating.
cmd_esl_seqstat = ["esl-seqstat", my_reference_fasta]
subprocess.call(cmd_esl_seqstat)

#Use the cmscan program to annotate RNAs represented in Rfam 
#in the S. clavuligerus genome

cmd_cmscan = ["cmscan", 
               "-Z",  "17.088172", "--tblout", "Sclav-genome.tblout",
               "--cut_ga", "--rfam", "--cpu", "18",
               "Rfam_2024-04.cm", my_reference_fasta]
# Open the output file for writing
with open("Sclav-genome.cmscan", "w") as output_file:
    # Redirect stdout to the output file
    subprocess.call(cmd_cmscan, stdout=output_file)

#Move the alignments to the "Infernal" directory
infernal_results = [infernal for infernal in os.listdir(os.getcwd()) if (infernal.endswith(".cmscan") or infernal.endswith(".tblout"))]  
for file_name in infernal_results:
    full_file_name = os.path.join(os.getcwd(), file_name)
    shutil.move(full_file_name, "results_sRNA/Infernal")

#Use the cmscan program to annotate RNAs represented in Rfam 
#in the conserved IGRs
cmd_esl_seqstat = ["esl-seqstat", IGRs_conserved_reference]
subprocess.call(cmd_esl_seqstat)
cmd_cmscan = ["cmscan", 
              "-Z",  "0.200716", "--tblout", "Sclav-conserved_IGRs.tblout",
              "--cut_ga", "--rfam", "--cpu", "18",
              "Rfam_2024-04.cm", IGRs_conserved_reference]
# Open the output file for writing
with open("Sclav-conserved_IGRs.cmscan", "w") as output_file:
    # Redirect stdout to the output file
    subprocess.call(cmd_cmscan, stdout=output_file)

#Move the alignments to the "Infernal" directory
infernal_results = [infernal for infernal in os.listdir(os.getcwd()) if (infernal.endswith(".cmscan") or infernal.endswith(".tblout"))]  
for file_name in infernal_results:
    full_file_name = os.path.join(os.getcwd(), file_name)
    shutil.move(full_file_name, "results_sRNA/Infernal")


#Use the cmscan program to annotate RNAs represented in Rfam 
#in the conserved IGRs with RNAz high class probability
cmd_esl_seqstat = ["esl-seqstat", my_conserved_IGRs_RNAz]
subprocess.call(cmd_esl_seqstat)
cmd_cmscan = ["cmscan", 
              "-Z",  "0.10411", "--tblout", "Sclav-conserved_IGRs_RNAz.tblout",
              "--cut_ga", "--rfam", "--cpu", "18",
              "Rfam_2024-04.cm", my_conserved_IGRs_RNAz]
               
# Open the output file for writing
with open("Sclav-conserved_IGRs_RNAz.cmscan", "w") as output_file:
    # Redirect stdout to the output file
    subprocess.call(cmd_cmscan, stdout=output_file)

#Move the alignments to the "Infernal" directory
infernal_results = [infernal for infernal in os.listdir(os.getcwd()) if (infernal.endswith(".cmscan") or infernal.endswith(".tblout"))]  
for file_name in infernal_results:
    full_file_name = os.path.join(os.getcwd(), file_name)
    shutil.move(full_file_name, "results_sRNA/Infernal")


infernal_annotations = pd.read_csv("results_sRNA/Infernal/Sclav-conserved_IGRs.tblout", comment='#', header=None)
infernal_annotations = infernal_annotations[0].str.split(expand=True)
infernal_annotations['description_of_target'] = infernal_annotations.apply(lambda row: ' '.join(str(val) for val in row[17:-1] if val is not None), axis=1)


infernal_annotations_filtered = infernal_annotations.iloc[:, [-1, 1,2,7,8,9,14,15]]


infernal_columns = ["target_name",  "accession", "query_name", "seq_from",
                    "seq_to", "strand", "score", "E-value"]
infernal_annotations_filtered.columns = infernal_columns

infernal_annotations_filtered.index.duplicated()

infernal_annotations_filtered = infernal_annotations_filtered.drop_duplicates(subset=['query_name'])



infernal_annotations_filtered = infernal_annotations_filtered.set_index('query_name')


putative_sRNAs_rfam_blast_infernal = pd.concat([putative_sRNAs_rfam_blast, 
                                                infernal_annotations_filtered], axis = 1)



putative_sRNAs_rfam_blast_infernal.to_csv ('results_sRNA/reference_results/putative_sRNAs_rfam_blast_infernal.csv', index = True, header = True)



end_time = time.time()
elapsed_time = end_time - start_time




#%% Define the coding potential of IGRs

#For this use CPC2 available at https://cpc2.gao-lab.org/index.php

#Call coding potential  results from CPC2
coding_potential = pd.read_table("results_sRNA/CPC2/conserved_IGRs_my_reference_cpc2.txt", sep="\t")
coding_potential_reordered = coding_potential[['#ID', 'label', 'peptide_length', 'Fickett_score', 'pI', 'ORF_integrity', 'coding_probability']]

coding_potential_reordered.set_index('#ID', inplace=True)


putative_sRNAs_rfam_blast_infernal_CP = pd.concat([putative_sRNAs_rfam_blast_infernal, 
                                                coding_potential_reordered], axis = 1)


putative_sRNAs_rfam_blast_infernal_CP.to_csv ('results_sRNA/reference_results/putative_sRNAs_rfam_blast_infernal_CP.csv', index = True, header = True)


#%% command line IntaRNA
"""
IntaRNA -t /media/usuario/lab_bioinformatica/Carlos_Caicedo-Montoya/sRNAs_Sclav/results_sRNA/reference_results/conserved_IGRs_my_reference_RNAz.dat.fasta -q /media/usuario/lab_bioinformatica/Carlos_Caicedo-Montoya/sRNAs_Sclav/raw_data/CDS_from_genomic_group1/GCF_005519465.1_ASM551946v1_cds_from_genomic.fasta --threads 6 --personality=IntaRNAsTar --out=IntaRNA_Sclav_RNAz_conserved IntaRNA -t /media/usuario/lab_bioinformatica/Carlos_Caicedo-Montoya/sRNAs_Sclav/results_sRNA/reference_results/conserved_IGRs_my_reference_RNAz.fasta -q /media/usuario/lab_bioinformatica/Carlos_Caicedo-Montoya/sRNAs_Sclav/raw_data/CDS_from_genomic_group1/GCF_005519465.1_ASM551946v1_cds_from_genomic.fasta --threads 6 --personality=IntaRNAsTar

IntaRNA -t /s_my_reference_RNAz.fasta -q /media/usuario/lab_bioinformatica/Carlos_Caicedo-Montoya/sRNAs_Sclav/raw_data/CDS_from_genomic_group1/GCF_005519465.1_ASM551946v1_cds_from_genomic.fasta --threads 6 --personality=IntaRNAsTar
"""


t_inta = "results_sRNA/reference_results/conserved_IGRs_my_reference_RNAz.fasta"
q_inta = "raw_data/CDS_from_genomic_group1/GCF_005519465.1_ASM551946v1_cds_from_genomic.fasta"


os.mkdir("results_sRNA/IntaRNA")
Inta_Results = os.path.join("results_sRNA/IntaRNA","IntaRNA_Sclav_RNAz_conserved")
cmd_inta = ["IntaRNA", '-t', t_inta, "-q", q_inta,  "--threads", "18", "--personality=IntaRNAsTar", "--out", Inta_Results]
subprocess.call(cmd_inta, text=True)


intarRNA_output = pd.read_table("results_sRNA/IntaRNA/IntaRNA_Sclav_RNAz_conserved", sep=";")

intarRNA_output_30 = intarRNA_output[intarRNA_output['E'] <= -30]

intarRNA_output_20 = intarRNA_output[intarRNA_output['E'] <= -20]

intarRNA_output_25 = intarRNA_output[intarRNA_output['E'] <= -25]

intarRNA_output_10 = intarRNA_output[intarRNA_output['E'] <= -10]



intarRNA_output_30.to_csv("results_sRNA/IntaRNA/IntaRNA_Sclav_RNAz_conserved_30", sep = "\t", header = True)
intarRNA_output_25.to_csv("results_sRNA/IntaRNA/IntaRNA_Sclav_RNAz_conserved_25", sep = "\t", header = True)
intarRNA_output_20.to_csv("results_sRNA/IntaRNA/IntaRNA_Sclav_RNAz_conserved_20", sep = "\t", header = True)
intarRNA_output_10.to_csv("results_sRNA/IntaRNA/IntaRNA_Sclav_RNAz_conserved_10", sep = "\t", header = True)




#Definir orientacion correcta4

#Course bash
#make the program executable

#df2--Promotech
#df3--Transtermhp
#1 +
#0 -


#%% General features of IGRS in S. clavuligerus



my_input_dir = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/raw_data/Genbank_group1"
my_reference = "GCF_005519465.1_ASM551946v1_genomic.gbff"

#length distribution call get_intergenic_regions without limiting the length of the IGRS
handle = os.path.join(my_input_dir, my_reference)

replicon_size = []
for i in SeqIO.parse(handle, "genbank"):
    replicon_size.append(len(i))
replicon_size = np.asarray(replicon_size)
genome_size = np.sum(replicon_size)


get_intergenic_regions(handle, intergenic_length=1, 
                       max_intergenic_length=10000000, 
                       output_dir="results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/reference_results")
my_reference_IGRS_fasta = os.path.join("results_sRNA_clade2_pangenome_gene_nontrimmed_50-600", 
                                       "reference_results", 
                                      my_reference.split(".")[0] + "_IGRs.fasta")
lengths_IGRs=[]
for i in SeqIO.parse(my_reference_IGRS_fasta, "fasta"):
    lengths_IGRs.append(len(i))
lengths_IGRs=np.asarray(lengths_IGRs)
lengths_IGRs_sorted = np.sort(lengths_IGRs)
min_IGR = np.min(lengths_IGRs)
max_IGR = np.max(lengths_IGRs)
median_lengths_IGRs= np.median(lengths_IGRs)
total_lengths_IGRs = np.sum(lengths_IGRs)

get_intergenic_regions(handle, intergenic_length=50,
                       max_intergenic_length=600,
                       output_dir="results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/reference_results")
my_reference_IGRS_fasta = os.path.join("results_sRNA_clade2_pangenome_gene_nontrimmed_50-600", 
                                       "reference_results", 
                                       my_reference.split(".")[0] + "_IGRs.fasta")
lengths_IGRs_reduced=[]
for i in SeqIO.parse(my_reference_IGRS_fasta, "fasta"):
    lengths_IGRs_reduced.append(len(i))

lengths_IGRs_reduced=np.asarray(lengths_IGRs_reduced)
lengths_IGRs_reduced_sorted = np.sort(lengths_IGRs_reduced)
min_IGR_reduced = np.min(lengths_IGRs_reduced)
max_IGR_reduced = np.max(lengths_IGRs_reduced)
median_lengths_IGRs_reduced= np.median(lengths_IGRs_reduced)

#Frequency of the type of IGRs 

my_reference_IGRS_fasta = os.path.join("results_sRNA_clade2_pangenome_gene_nontrimmed_50-600", 
                                       "IGRs", 
                                       my_reference.split(".")[0] + "_IGRs_reduced.fasta")

flanking_genes = []
for i in SeqIO.parse(my_reference_IGRS_fasta, "fasta"):
    flanking_genes.append(i.id.split("++")[3])


# using Counter to find frequency of elements
frequency = collections.Counter(flanking_genes)

# Convert frequency to a dictionary an then to a dataframe
frequency = dict(frequency)
df = pd.DataFrame.from_dict(frequency, orient='index')
df.rename(columns={0: "flanking_genes"}, inplace=True)
labels = list(df.index)

#determine GC content of intergenic regions

GCs = []
file = "results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/IGRs/GCF_005519465_IGRs_reduced.fasta"
for seq_record in SeqIO.parse(file, "fasta"):
    seq_record.seq
    GC_perc =  gc_fraction(seq_record.seq)
    GCs.append(GC_perc)


GCs = np.asarray(GCs)
GCs_sorted = np.sort(GCs)
min_GCs = np.min(GCs)
max_GCs = np.max(GCs)
median_GCs= np.median(GCs)

GC_mean = np.mean(GCs)


f, ((ax1, ax2), (ax3, ax4))= plt.subplots(2, 2)
#plt.f(figsize=(8, 6))
ax1.hist(np.asarray(lengths_IGRs), bins=50, color = "darkcyan")
ax1.grid()
ax1.set_xlabel("Length of IGR (bp)")
ax1.set_ylabel('Number of IGRs', size = 12)
values1 = np.array([0,0.25,0.5,0.75,1])
x1 = np.quantile(np.asarray(lengths_IGRs), values1)

ax2.hist(np.asarray(lengths_IGRs_reduced), color  ="Deepskyblue")
ax2.grid()
ax2.set_xlabel("Length of IGR (bp)")
values2 = np.array([0,0.25,0.5,0.75,1])
x2 = np.quantile(np.asarray(lengths_IGRs_reduced), values2)

#Plot the frequencies as barplots
values3 = df.iloc[:, 0]
values3 = np.asarray(values3)
x3 = np.arange(len(labels))  # the label locations
ax3.bar(x3, values3, tick_label = labels, color = "steelblue")
ax3.set_ylabel('Number of IGRs', size = 12)
ax3.grid()

ax4.hist(np.asarray(GCs), color = "mediumslateblue")
ax4.grid()
ax4.set_xlabel("% GC")
values4 = np.array([0,0.25,0.5,0.75,1])
x4 = np.quantile(np.asarray(lengths_IGRs), values1)


f.savefig ("results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/reference_results/F1_A-B-C-D.pdf",
           dpi=300, format = "pdf")

#%% count sequences in cluster of conserved IGRs fasta file

# Path to the folder where the fasta files are located
folder_path = "results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/IGRs_cluster"

sequences_count_per_file = [(file, sum(1 for _ in SeqIO.parse(os.path.join(folder_path, file), 'fasta')))
                            for file in os.listdir(folder_path) if file.endswith('.fasta')]


# Print the number of sequences per file
sequences_per_cluster=[]
indentifiers_per_cluster = []
for file, num_sequences in sequences_count_per_file:
    sequences_per_cluster.append(num_sequences)
    indentifiers_per_cluster.append(file)


sequences_per_cluster = np.asarray(sequences_per_cluster)
f2, ax1= plt.subplots()
ax1.hist(np.asarray(sequences_per_cluster), bins = 20, color = "indigo")
ax1.grid()
ax1.set_xlabel("Number of sequences", size = 12)
ax1.set_ylabel('Frequency', size = 12)
values1 = np.array([0,0.25,0.5,0.75,1])
x1 = np.quantile(np.asarray(sequences_per_cluster), values1)

unique_values, frequencies = np.unique(sequences_per_cluster, return_counts=True)

frequencies_conserved_IGRs = dict(zip(unique_values, frequencies))


putative_sRNAs_rfam_blast_infernal_CP = pd.read_table("results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/reference_results/putative_sRNAs_rfam_blast_infernal_CP.csv", 
                                                      sep=",", index_col=0 )

sequences_count_per_file_df = pd.DataFrame(sequences_count_per_file, columns=['index', 'number_of_genomes'])

# Establecer la columna 'index' como el índice
sequences_count_per_file_df.set_index('index', inplace=True)
sequences_count_per_file_df.index = sequences_count_per_file_df.index.str.replace('.fasta', '')




putative_sRNAs_rfam_blast_infernal_CP_NOG = putative_sRNAs_rfam_blast_infernal_CP.join(sequences_count_per_file_df) 
                                                       
#Number of conserved IGRs with secondary structure conservation
RNAz_count = putative_sRNAs_rfam_blast_infernal_CP_NOG['RNAz'].sum()

#Number of conserved IGRs with a promoter
Promoter_count = putative_sRNAs_rfam_blast_infernal_CP_NOG['Promoter1'].count()
Promoter_count2 = putative_sRNAs_rfam_blast_infernal_CP_NOG['Promoter2'].count()

#Number of conserved IGRs with a promoter
Terminator_count = putative_sRNAs_rfam_blast_infernal_CP_NOG['Terminator1'].count()

#Number of conserved IGRs annotation in RFAM
annotation_count = putative_sRNAs_rfam_blast_infernal_CP_NOG['target_name'].count()

#number of conserved IGRs with coding potential

index_coding = putative_sRNAs_rfam_blast_infernal_CP_NOG.loc[putative_sRNAs_rfam_blast_infernal_CP_NOG['label'] == 'coding'].index[0]


#%% Extract RNAz  RNA-class probability


def extract_rnaz_info(file_path):
    # Define an empty list to store the extracted information
    results_list = []

    # Open the file and read it line by line
    with open(file_path, 'r') as file:
        lines = file.readlines()
        # Initialize variables to store the current reading direction and SVM RNA-class probability
        reading_direction = None
        svm_probability = None
        for line in lines:
            # Check if the line marks the beginning of a new RNAz 2.1 result section
            if line.startswith("############################  RNAz 2.1  ##############################"):
                # If the current reading direction and SVM RNA-class probability are not None, add them to the list
                if reading_direction is not None and svm_probability is not None:
                    results_list.append((reading_direction, svm_probability))
                # Reset the variables for the new result section
                reading_direction = None
                svm_probability = None
            # Extract the reading direction line
            elif "Reading direction:" in line:
                reading_direction = line.strip()
            # Extract the SVM RNA-class probability line
            elif "SVM RNA-class probability:" in line:
                svm_probability = line.strip()

    # Add the last extracted reading direction and SVM RNA-class probability to the list
    if reading_direction is not None and svm_probability is not None:
        results_list.append((reading_direction, svm_probability))

    return results_list


rnaz_results = [rnaz for rnaz in os.listdir("results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/RNAz_out") if rnaz.endswith(".out") ]
rnaz_results.remove("RNAz_all.out")

folder_rnaz = os.path.join(os.getcwd(),"results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/RNAz_out")
results_dict_rnaz = {}

# Iterate over the list of files
for file_name in rnaz_results:
    # Get the file name without the ".fasta.out" extension
    key = file_name.replace('.fasta.out', '')
    file_path = os.path.join(folder_rnaz, file_name)
    values = extract_rnaz_info(file_path)
    # Add to the dictionary
    results_dict_rnaz[key] = values

results_dict_rnaz = {key: value for key, value in results_dict_rnaz.items() if value}

results_dict_rnaz_df = pd.DataFrame.from_dict(results_dict_rnaz, orient='index', columns=['Strand1', 'Strand2'])

# Split the columns containing tuples into two separate columns
results_dict_rnaz_df[['Strand1', 'RNA_Class_probability1']] = results_dict_rnaz_df['Strand1'].apply(pd.Series)
results_dict_rnaz_df[['Strand2', 'RNA_Class_probability2']] = results_dict_rnaz_df['Strand2'].apply(pd.Series)
results_dict_rnaz_df = results_dict_rnaz_df[['Strand1', 'RNA_Class_probability1', 'Strand2', 'RNA_Class_probability2']]
results_dict_rnaz_df['Strand1'] = results_dict_rnaz_df['Strand1'].str.replace('Reading direction: ', '')
results_dict_rnaz_df['Strand2'] = results_dict_rnaz_df['Strand2'].str.replace('Reading direction: ', '')
results_dict_rnaz_df['RNA_Class_probability1'] = results_dict_rnaz_df['RNA_Class_probability1'].str.replace('SVM RNA-class probability: ', '')
results_dict_rnaz_df['RNA_Class_probability2'] = results_dict_rnaz_df['RNA_Class_probability2'].str.replace('SVM RNA-class probability: ', '')

position = putative_sRNAs_rfam_blast_infernal_CP_NOG.columns.get_loc('RNAz') + 1  # Obtener la posición después de 'RNAz'
for col_name, col_data in results_dict_rnaz_df.items():
    putative_sRNAs_rfam_blast_infernal_CP_NOG.insert(position, col_name, col_data)

putative_sRNAs_rfam_blast_infernal_CP.to_csv ('results_sRNA/reference_results/putative_sRNAs_rfam_blast_infernal_CP.csv', index = True, header = True)



#%% Analyze in detail the IGRs

putative_sRNAs_summary = putative_sRNAs_rfam_blast_infernal_CP_NOG[[ 'Promoter1', 
                                                                    'Terminator1',
                                                                    'RNAz',
                                                                    'Strand1',
                                                                    'RNA_Class_probability1',
                                                                    'Strand2',
                                                                    'RNA_Class_probability2',
                                                                    'subject',
                                                                    'target_name',
                                                                    'number_of_genomes', 
                                                                    'label']].copy()


#List for Venn Diagram
index_rnaz = putative_sRNAs_summary.loc[putative_sRNAs_summary['RNAz'] == True].index
list_index_RNAz = index_rnaz.tolist()


index_promoter = putative_sRNAs_summary['Promoter1'].dropna().index
list_index_promoter = index_promoter.tolist()


index_terminator = putative_sRNAs_summary['Terminator1'].dropna().index
list_index_terminator = index_terminator.tolist()


has_infernal = putative_sRNAs_summary['target_name'].dropna().index
list_has_infernal = has_infernal.tolist()

#consider if previous annotations of sRNAS were detected in your analysis

file_gff_all = os.path.join(os.getcwd(), "raw_data", 
                            "S-clavuligerus_ATCC27064_annotations", 
                            "GCF_005519465.1_ASM551946v1_genomic.gff")


gff_all = pd.read_csv(file_gff_all, comment='#', 
                            sep = "\t", header=None)

gff_columns = ["seqname", "source", "feature", "start", "end", "score", 
                   "strand", "frame", "attribute" ]
gff_all.columns = gff_columns
gff_all_cmsearch = gff_all[gff_all['source'] == 'cmsearch']
del gff_all

#Read Infernal annotations for all the genome
infernal_annotations = pd.read_csv("results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/Infernal/Sclav-genome.tblout", comment='#', header=None, delimiter='\t')
infernal_annotations = infernal_annotations[0].str.split(expand=True)
infernal_annotations['description_of_target'] = infernal_annotations.apply(lambda row: ' '.join(str(val) for val in row[17:-1] if val is not None), axis=1)

infernal_annotations_filtered = infernal_annotations.iloc[:, [-1, 1,2,7,8,9,14,15]]
infernal_columns = ["target_name",  "accession", "query_name", "seq_from",
                    "seq_to", "strand", "score", "E-value"]
infernal_annotations_filtered.columns = infernal_columns

infernal_annotations_filtered = infernal_annotations_filtered.drop_duplicates(subset=['seq_from'])
infernal_annotations_filtered = infernal_annotations_filtered.set_index('target_name')

infernal_annotations_filtered_sRNa = infernal_annotations_filtered.copy() 
infernal_annotations_filtered_sRNa = infernal_annotations_filtered_sRNa.drop('Bacterial large subunit ribosomal RNA')
infernal_annotations_filtered_sRNa = infernal_annotations_filtered_sRNa.drop('Archaeal large subunit ribosomal RNA')
infernal_annotations_filtered_sRNa = infernal_annotations_filtered_sRNa.drop('Bacterial small subunit ribosomal RNA')
infernal_annotations_filtered_sRNa = infernal_annotations_filtered_sRNa.drop('Eukaryotic large subunit ribosomal RNA')
infernal_annotations_filtered_sRNa = infernal_annotations_filtered_sRNa.drop('Archaeal small subunit ribosomal RNA')
infernal_annotations_filtered_sRNa = infernal_annotations_filtered_sRNa.drop('5S ribosomal RNA')
infernal_annotations_filtered_sRNa = infernal_annotations_filtered_sRNa.drop('tRNA')
infernal_annotations_filtered_sRNa = infernal_annotations_filtered_sRNa.drop('transfer-messenger RNA')
infernal_annotations_filtered_sRNa = infernal_annotations_filtered_sRNa.drop('Bacterial small signal recognition particle')
infernal_annotations_filtered_sRNa = infernal_annotations_filtered_sRNa.drop('Bacterial RNase P class A')
infernal_annotations_filtered_sRNa = infernal_annotations_filtered_sRNa.drop('Group II catalytic intron')






# See the distribution of flanking genes in conserved igrs

index_conserved_IGRs = putative_sRNAs_summary.index.tolist()
index_series = pd.Series(index_conserved_IGRs)

# List to store the results
Class_Conserved_IGRs = []

# Iterate over each element in the list
for string in index_conserved_IGRs:
    parts = string.split("++")  # Split the string by '++'
    last_part = parts[-1]  # Get the last part
    Class_Conserved_IGRs.append(last_part)  # Append to the new list

# Create the histogram
f3, ax1 = plt.subplots()
ax1.hist(Class_Conserved_IGRs, color = "navy")
ax1.grid()
ax1.set_ylabel('Frequency', size = 12)


frequency_class_conserved_IGRs = dict(Counter(Class_Conserved_IGRs))
frequency_class_conserved_IGRs = list(frequency_class_conserved_IGRs.values())
expected_freq = [0.25*len(Class_Conserved_IGRs), 
                 0.25*len(Class_Conserved_IGRs),
                 0.25*len(Class_Conserved_IGRs), 
                 0.25*len(Class_Conserved_IGRs)]
chi2, p_valor, _, _ = chi2_contingency([frequency_class_conserved_IGRs,
                                        expected_freq])

### MOTIFS in conserved IGRs

def format_promoters_in_IGRs(promoters):
    # Drop rows where 'Promoter1' is NaN
    promoters = promoters.dropna(subset=['Promoter1'])
    # Ensure 'Promoter1' is a string before splitting
    promoters.loc[:, 'Promoter1'] = promoters['Promoter1'].apply(lambda x: ','.join(map(str, x)))
    # Split 'Promoter1' into separate columns
    promoters[['start', 'end', 'strand', 'sequence']] = promoters['Promoter1'].str.split(',', expand=True)   
    # Drop columns starting with 'Promoter'
    columns_to_drop = promoters.filter(regex='^Promoter').columns
    promoters = promoters.drop(columns=columns_to_drop)
    
    return promoters

promoters_in_IGRs_conserved = format_promoters_in_IGRs(promoters_in_IGRs_conserved)


### MOTIF in promoter in all IGRs
my_assembly_name1= my_assembly_name(my_reference)
fasta_file_IGRs = f"results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/IGRs/{my_assembly_name1}_IGRs_reduced.fasta"
output_file_name = "Promoters_in_ALL_IGRs.csv"  
promoters_in_ALL_IGRs = promoters_in_IGRs(fasta_file_IGRs, promoters_promotech, output_file_name)
promoters_in_ALL_IGRs = format_promoters_in_IGRs(promoters_in_ALL_IGRs)




def write_fasta_file(df, fasta_file):
    with open(fasta_file, 'w') as f:
        for index, row in df.iterrows():
            fasta_id = index  # Use the index of each row as the fasta_id
            sequence = row['sequence']
            f.write(f'>{fasta_id}\n{sequence}\n')
            
write_fasta_file(promoters_in_IGRs_conserved, 
                 'results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/reference_results/promoters_in_IGRs_conserved.fasta')
write_fasta_file(promoters_in_ALL_IGRs, 
                 'results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/reference_results/promoters_in_all_IGRs.fasta')




#trasntermHP for all IGRs

terminator_results_ALL_IGRS_path = os.path.join("results_sRNA", 
                                                      "TransTermHP", 
                                                      "terminator_in_All_IGRs")
terminator_in_all_IGRs = terminator_in_IGRs(fasta_file_IGRs, 
                                        terminators_transtermhp,                                         
                                        terminator_results_ALL_IGRS_path)


all_IGRs_term_prom = pd.concat([promoters_in_ALL_IGRs, 
                                                terminator_in_all_IGRs], axis = 1)

all_IGRs_term_prom.to_csv ('results_sRNA/reference_results/all_IGRs_term_prom.csv', index = True, header = True)

#BLASTn search
my_IGRs_reference = "results_sRNA_clade2_pangenome_gene_nontrimmed_50/IGRs/GCF_005519465_IGRs_reduced.fasta"
rfam_db = os.path.join(os.getcwd(), "raw_data", "rfam_sequences", "Rfam.fa")


cline_rfam = NcbiblastnCommandline(query=my_IGRs_reference, 
                                  db=rfam_db,
                                  evalue=1e-6, 
                                  out = os.path.join("results_sRNA", "reference_results", "my_conserved_IGRs_RNAz_blastn_out.txt"), 
                                  task = 'blastn',
                                  num_threads = 10,
                                  max_hsps = 1,
                                  max_target_seqs=1,
                                  #dust = 'no',
                                  outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue")
stdout, stderr = cline_rfam()


rfam_results = pd.read_table(os.path.join("results_sRNA", "reference_results", "my_conserved_IGRs_RNAz_blastn_out.txt"), sep='\t', header=None)

rfam_results.columns = ["query", "subject", "identity", "coverage",
                   "qlength", "slength", "alength",
                   "bitscore", "E-value"]
rfam_results.set_index("query", inplace=True)



#Use the cmscan program to annotate RNAs represented in Rfam 
#in all  IGRs
cmd_esl_seqstat = ["esl-seqstat", fasta_file_IGRs]
subprocess.call(cmd_esl_seqstat)
cmd_cmscan = ["cmscan", 
              "-Z",  "1.924788", "--tblout", "Sclav-all_IGRs.tblout",
              "--cut_ga", "--rfam", "--cpu", "18",
              "Rfam_2024-04.cm", fasta_file_IGRs]
# Open the output file for writing
with open("Sclav-all_IGRs.cmscan", "w") as output_file:
    # Redirect stdout to the output file
    subprocess.call(cmd_cmscan, stdout=output_file)

#Move the alignments to the "Infernal" directory
infernal_results = [infernal for infernal in os.listdir(os.getcwd()) if (infernal.endswith(".cmscan") or infernal.endswith(".tblout"))]  
for file_name in infernal_results:
    full_file_name = os.path.join(os.getcwd(), file_name)
    shutil.move(full_file_name, "results_sRNA/Infernal")


infernal_annotations = pd.read_csv("results_sRNA/Infernal/Sclav-all_IGRs.tblout", comment='#', header=None)
infernal_annotations = infernal_annotations[0].str.split(expand=True)
infernal_annotations['description_of_target'] = infernal_annotations.apply(lambda row: ' '.join(str(val) for val in row[17:-1] if val is not None), axis=1)


infernal_annotations_filtered = infernal_annotations.iloc[:, [-1, 1,2,7,8,9,14,15]]


infernal_columns = ["target_name",  "accession", "query_name", "seq_from",
                    "seq_to", "strand", "score", "E-value"]
infernal_annotations_filtered.columns = infernal_columns
infernal_annotations_filtered.index.duplicated()
infernal_annotations_filtered = infernal_annotations_filtered.drop_duplicates(subset=['query_name'])
infernal_annotations_filtered = infernal_annotations_filtered.set_index('query_name')
all_IGRs_infernal_prom_term = pd.concat([all_IGRs_term_prom, 
                                                infernal_annotations_filtered], axis = 1)

all_IGRs_infernal_prom_term.to_csv ('results_sRNA/reference_results/all_IGRs_infernal_prom_term.csv', index = True, header = True)

all_IGRs_infernal_prom_term.drop(columns=["Terminator2",
                                          "Terminator3", 
                                          "Terminator4", 
                                          "Terminator5", 
                                          "Terminator6",
                                          "Terminator7"], inplace = True)


#number of promoters in all IGRs

number_promoters_all_IGRs = all_IGRs_infernal_prom_term['sequence'].notna().sum()

number_terminators_all_IGRs = all_IGRs_infernal_prom_term['Terminator1'].notna().sum()








#%% length distribution of alignments of conserved IGRS


msa_conserved_IGRs = [msa for msa in os.listdir("results_sRNA/IGRs_alignments") if msa.endswith(".aln") ]

lengths_msa_conserved_IGRs=[]
for msa_files in msa_conserved_IGRs:
    msa = AlignIO.read(os.path.join("results_sRNA", "IGRs_alignments",  msa_files), "clustal")
    lengths_msa_conserved_IGRs.append(msa.get_alignment_length())
               
f1, ax= plt.subplots(1, 2)
ax[0].boxplot(np.asarray(lengths_msa_conserved_IGRs))
ax[0].grid()
ax[0].set_ylabel("Length of MSA (bp)")

ax[1].hist(np.asarray(lengths_msa_conserved_IGRs))
ax[1].grid()
ax[1].set_xlabel("Length of MSA (bp)")
ax[1].set_ylabel("Number of MSA")
values = np.array([0,0.25,0.5,0.75,1])
x = np.quantile(np.asarray(lengths_msa_conserved_IGRs), values)














 
import subprocess 
i = "conserved_IGRs_my_reference.fasta"
num = len([1 for line in open(i) if line.startswith(">")])    
print(num, "IGRs in file", i)   

s = [aln for aln in os.listdir("results_sRNA/IGRs_cluster") if aln.endswith(".aln") ]  


#Check the output of transtermhp and promoteh by opening  the region that contains 
#the terminator or the promoter

record_iterator = SeqIO.parse(my_reference_fasta, "fasta")
first_record = next(record_iterator)
terminator_example = first_record.seq[63741:63775] 
terminator_example.reverse_complement()

reverse_complement = (rec.reverse_complement(id=rec.id, description = "reverse_complement")
           for rec in SeqIO.parse(my_reference_fasta, "fasta"))
SeqIO.write(reverse_complement, "results_sRNA/reference_results/" + my_assembly_name(my_reference) + "_rev_comp.fasta", "fasta")

 






