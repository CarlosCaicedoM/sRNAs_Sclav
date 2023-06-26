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
#pip install wget
# conda install -c bioconda clustalw
# conda install -c bioconda rnaz
# conda install -c bioconda infernal

#%%  Functions
from IPython import get_ipython
get_ipython().magic('reset -sf')
#run_line_magic(magic_name, parameter_s).get_ipython().magic('reset -sf')


import sys
import Bio
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
import re
import wget
import gzip

### Define function to extract IGRs
def get_intergenic_regions(handle, intergenic_length, max_intergenic_length, output_dir = "results_sRNA/IGRs"):   

    intergenic_records = []
    
    for seq_record in SeqIO.parse(handle, "genbank"):
        gene_locations = []
        for feature in seq_record.features:
            if feature.type == 'gene':
                my_start = feature.location.start.position
                my_end = feature.location.end.position
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
    my_new_reference  = my_input_dir + my_reference
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
        my_file = output_dir + "/" + i 
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
        
        clustal_cline = ClustalwCommandline("clustalw2", infile=output_file, 
                                         outfile=output_file.split()[0] + ".aln", 
                                         outorder="INPUT")
        stdout, stderr = clustal_cline()
        
        
        remove(id_file)
        

def my_assembly_name(my_reference):
    my_new_reference  = my_input_dir + my_reference
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

from Bio import SeqIO
import os
import sys
from random import randint

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

#Extract IGRs from the genomes to compare

#Define a folder to store temporary files and intermediate results


my_current_dirs = [name for name in os.listdir(".") if os.path.isdir(name)]




if "results_sRNA" in my_current_dirs:
    shutil.rmtree("results_sRNA")
    os.mkdir("results_sRNA")
else:
    os.mkdir("results_sRNA")

"""
for dir_path in my_current_dirs:
    try:
        shutil.rmtree(dir_path)
    except OSError as e:
        print("Error: %s : %s" % (dir_path, e.strerror))
"""


os.mkdir(os.path.join("results_sRNA", "IGRs_blast"))
os.mkdir(os.path.join("results_sRNA", "IGRs_cluster"))
os.mkdir(os.path.join("results_sRNA", "IGRs_alignments"))
os.mkdir(os.path.join("results_sRNA", "reference_results"))  
os.mkdir(os.path.join("results_sRNA", "IGRs"))


#Define the folder with the genomes to analyze
#my_input_dir = "/media/usuario/lab_bioinformatica/Carlos_Caicedo-Montoya/sRNAs_Sclav/raw_data/Genbank_group1/"
my_input_dir = "raw_data/Genbank_group1/"
#my_input_dir = "/media/usuario/lab_bioinformatica/Carlos_Caicedo-Montoya/sRNAs_Sclav/test_Actinoplanes/"


#to get the data consider to use the following comand

"""
conda install -c bioconda ncbi-genome-download
ncbi-genome-download bacteria -F genbank -A Streptomyces_Assemby_Accession_Group1_Sclav.txt -v

#To unzip multiple files
chmod +x unzip_several.sh
./unzip_several.sh

"""


#Extract IGRs for the genomes in my input directory
my_files = os.listdir(my_input_dir)
for i in my_files:
    print("Processing file %s" %i )
    get_intergenic_regions(os.path.join(my_input_dir, i), intergenic_length=130, max_intergenic_length=600)

#Summarize the number of IGRs in each genome
my_results = os.listdir(os.path.join("results_sRNA", "IGRs"))
for i in my_results:
    num = len([1 for line in open(os.path.join("results_sRNA", "IGRs", i)) if line.startswith(">")])    
    print(num, "IGRs in file", i)                                

#Cut the borders of the sequence to avoid to include regulatory elements in the analysis
#that can belong to the 5' or 3' UTRs
for i in my_results:
   print("Processing file %s" %i )
   trim_sequences(os.path.join("results_sRNA", "IGRs", i), 40, 40)
   remove(os.path.join("results_sRNA", "IGRs", i))


## Create the databases

#my_reference = "GCF_001553785.1_ASM155378v1_genomic.gbff"  
my_reference = "GCF_005519465.1_ASM551946v1_genomic.gbff"

#install blast + suite
make_blast_db()

         
### Run reciprocal best blast 
threads = 8
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

#muscle_exe = "./muscle"
#makemuscle executable
#subprocess.call(["chmod", "+x", "muscle"])
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

headers_IGRs_conserved = []
flanking_genes = []
coordinates = []
orientation  = []
for seq_record in SeqIO.parse(IGRs_conserved_reference, "fasta"):
    headers_IGRs_conserved.append(seq_record.id)
    flanking_genes.append(seq_record.id.split("++")[1])
    coordinates.append(seq_record.id.split("++")[2])
    orientation.append(seq_record.id.split("++")[3])

df1 = pd.DataFrame({'headers':headers_IGRs_conserved, 
                    'flanking_genes':flanking_genes,
                    'coordinates':coordinates, 
                    'orientation':orientation})

IGRs_conserved_dict = {'headers':headers_IGRs_conserved, 
                    'flanking_genes':flanking_genes,
                    'coordinates':coordinates, 
                    'orientation':orientation}

IGRs_conserved_dict2 = dict(zip(headers_IGRs_conserved, coordinates))


#Locate promoters predicted by G4PromFinder or Promotech in the conserved IGRs
empty_lists =  [[] for _ in range(len(coordinates))]
promoters_in_IGRs = {headers_IGRs_conserved[i]:empty_lists[i] for i in range(len(headers_IGRs_conserved))}

print("Locating promoter in IGRs sequences")

for key, value in IGRs_conserved_dict2.items():
    start_IGR = int(IGRs_conserved_dict2[key].split("-")[0])
    end_IGR = int(IGRs_conserved_dict2[key].split("-")[1])
    for j in promoters_promotech.itertuples():
        if start_IGR <= j.start and end_IGR >= j.end and key.split("++")[0] == j.chrom:
            promoter_region = (j.start, j.end, j.strand)
            promoters_in_IGRs[key].append(promoter_region)
    

without_promoter = []
for key, value in promoters_in_IGRs.items():
    if len(value)==0:
        without_promoter.append(1)
        
number_of_promoters = []
for key, value in promoters_in_IGRs.items():
    number_of_promoters.append(len(value))


#Save promoters in IGRs
df2= pd.DataFrame.from_dict(promoters_in_IGRs, orient="index")
df2.to_csv ('results_sRNA/reference_results/Promoters_in_conserved_IGRs.csv', index = True, header = False)


#Determine terminators in the genome
os.mkdir(os.path.join("results_sRNA", "TransTermHP"))

#This function was inspired by the code found here (although modified for the purposes)
#of the current approach
#https://github.com/galaxyproject/tools-iuc/blob/master/tools/transtermhp/transtermhp.py


#convert genbank format to ptt
file = os.path.join(my_input_dir + my_reference)
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
df3.to_csv ('results_sRNA/reference_results/Terminators_in_conserved_IGRs.csv', index = True, header = False)

#Check the number of conserved IGRs with terminator 
has_terminators= []
for key, value in terminators_in_IGRs.items():
    if len(value) != 0:
        has_terminators.append(value)

without_terminator = []
for key, value in terminators_in_IGRs.items():
    if len(value)==0:
        without_terminator.append(1)

with_terminator = len(terminators_in_IGRs) - len(without_terminator)


#Run RNAz using the conserved IGR sequences
subprocess.call('./run_RNAz_local.sh')
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

#Download rfam database
ulr_rfam_fasta = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/Rfam.fa.gz"
filename_rfam_fasta = wget.download(ulr_rfam_fasta)


with gzip.open(filename_rfam_fasta, 'rb') as f_in:
    with open('Rfam.fa', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

#Move the downloaded Rfam database to the corresponding folder
os.mkdir(os.path.join("raw_data", "rfam_sequences"))
shutil.move("Rfam.fa", os.path.join(os.getcwd(), "raw_data", "rfam_sequences"))


#make Rfam.fa as blast db
cline_rfam_db = NcbimakeblastdbCommandline(dbtype="nucl", 
                                   input_file=rfam_db, 
                                   out = rfam_db)
stdout, stderr = cline_rfam_db()


#BLASTn search
cline_rfam = NcbiblastnCommandline(query=my_conserved_IGRs_RNAz, 
                                  db=rfam_db,
                                  evalue=1e-6, 
                                  out = os.path.join("results_sRNA", "reference_results", "my_conserved_IGRs_RNAz_blastn_out.txt"), 
                                  task = 'blastn',
                                  num_threads = 1,
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

test = pd.concat([df2, df3], axis = 1)


rnaz_ids = []
for rnaz in records_rnaz:
    rnaz_ids.append(rnaz.id)

have_rnaz =  [True for _ in range(len(rnaz_ids))]

rnaz_dict = dict(zip(rnaz_ids, have_rnaz))
rnaz_df = pd.DataFrame.from_dict(rnaz_dict, orient="index")

test2 = pd.concat([test, rnaz_df], axis = 1)


putative_sRNAs_rfam_blast = pd.concat([test2, rfam_results], axis = 1)


#Anotar con RFAM
url = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz'
filename = wget.download(url)
url2="ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin"
filename2  = wget.download(url2)


with gzip.open(filename, 'rb') as f_in:
    with open('Rfam.cm', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


os.mkdir("results_sRNA/Infernal")
Infernal_Results = os.path.join("results_sRNA/Infernal","Infernal_Sclav")
cmd_cmpress = ["cmpress", 'Rfam.cm']
subprocess.call(cmd_cmpress)

#corregir gunzip con python

esl-seqstat infernal-1.1.2/mrum-genome.fa

cmscan -Z 5.874406 --cut_ga --rfam --nohmmonly --tblout mrum-genome.tblout --fmt 2 --clanin Rfam.clanin Rfam.cm tutorial/mrum-genome.fa > mrum-genome.cmscan


#Definir orientacion correcta
#Anotar conserved IGRs not just the ones with RNAz results
#Promotech y trasntermHP para todos los IGRs


#df2--Promotech
#df3--Transtermhp
#1 +
#0 -



#%%Pending tasks
##command line IntaRNA
"""
IntaRNA -t /media/usuario/lab_bioinformatica/Carlos_Caicedo-Montoya/sRNAs_Sclav/results_sRNA/reference_results/conserved_IGRs_my_reference_RNAz.dat.fasta -q /media/usuario/lab_bioinformatica/Carlos_Caicedo-Montoya/sRNAs_Sclav/raw_data/CDS_from_genomic_group1/GCF_005519465.1_ASM551946v1_cds_from_genomic.fasta --threads 6 --personality=IntaRNAsTar --out=IntaRNA_Sclav_RNAz_conserved IntaRNA -t /media/usuario/lab_bioinformatica/Carlos_Caicedo-Montoya/sRNAs_Sclav/results_sRNA/reference_results/conserved_IGRs_my_reference_RNAz.fasta -q /media/usuario/lab_bioinformatica/Carlos_Caicedo-Montoya/sRNAs_Sclav/raw_data/CDS_from_genomic_group1/GCF_005519465.1_ASM551946v1_cds_from_genomic.fasta --threads 6 --personality=IntaRNAsTar

IntaRNA -t /s_my_reference_RNAz.fasta -q /media/usuario/lab_bioinformatica/Carlos_Caicedo-Montoya/sRNAs_Sclav/raw_data/CDS_from_genomic_group1/GCF_005519465.1_ASM551946v1_cds_from_genomic.fasta --threads 6 --personality=IntaRNAsTar
"""


t_inta = "results_sRNA/reference_results/conserved_IGRs_my_reference_RNAz.fasta"
q_inta = "raw_data/CDS_from_genomic_group1/GCF_005519465.1_ASM551946v1_cds_from_genomic.fasta"


os.mkdir("results_sRNA/IntaRNA")
Inta_Results = os.path.join("results_sRNA/IntaRNA","IntaRNA_Sclav_RNAz_conserved")
cmd_inta = ["IntaRNA", '-t', t_inta, "-q", q_inta,  "--threads", "6", "--personality=IntaRNAsTar", "--out", Inta_Results]
subprocess.call(cmd_inta, text=True)


intarRNA_output = pd.read_table("results_sRNA/IntaRNA/IntaRNA_Sclav_RNAz_conserved", sep=";")

intarRNA_output_30 = intarRNA_output[intarRNA_output['E'] <= -30]

intarRNA_output_20 = intarRNA_output[intarRNA_output['E'] <= -20]

intarRNA_output_25 = intarRNA_output[intarRNA_output['E'] <= -25]

intarRNA_output_10 = intarRNA_output[intarRNA_output['E'] <= -10]



intarRNA_output_30.to_csv("IntaRNA_Sclav_RNAz_conserved_30", sep = "\t", header = True)
intarRNA_output_25.to_csv("IntaRNA_Sclav_RNAz_conserved_25", sep = "\t", header = True)
intarRNA_output_20.to_csv("IntaRNA_Sclav_RNAz_conserved_20", sep = "\t", header = True)
intarRNA_output_10.to_csv("IntaRNA_Sclav_RNAz_conserved_10", sep = "\t", header = True)





 #Correct rbbh() directory for the output 

#Time of the program to run

#Delete blast databases

#Course bash

#make the program executable
#%% General features of IGRS in S. clavuligerus

#length distribution of conserved IGRS
lengths_conserved_IGRs=[]
for i in SeqIO.parse(IGRs_conserved_reference, "fasta"):
    lengths_conserved_IGRs.append(len(i))


f, ax= plt.subplots(1, 2)
ax[0].boxplot(np.asarray(lengths_conserved_IGRs))
ax[0].grid()
ax[0].set_ylabel("Length of IGR (bp)")

ax[1].hist(np.asarray(lengths_conserved_IGRs))
ax[1].grid()
ax[1].set_xlabel("Length of IGR (bp)")
ax[1].set_ylabel("Number of IGRs")
values = np.array([0,0.25,0.5,0.75,1])
x = np.quantile(np.asarray(lengths_conserved_IGRs), values)

#f.savefig ("reference_results/IGRs_boxplot.pdf", dpi=1200, format = "pdf")

#length distribution of alignments of conserved IGRS


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



#Frequency of the type of IGRs 

my_reference_IGRS_fasta = os.path.join("IGRs", my_reference.split(".")[0] + "_IGRs_reduced.fasta")

flanking_genes = []
for i in SeqIO.parse(my_reference_IGRS_fasta, "fasta"):
    flanking_genes.append(i.id.split("++")[3])


# importing the module
import collections
import numpy as np

# using Counter to find frequency of elements
frequency = collections.Counter(flanking_genes)

# Convert frequency to a dictionary an then to a dataframe
frequency = dict(frequency)
df = pd.DataFrame.from_dict(frequency, orient='index')
df.rename(columns={0: "flanking_genes"}, inplace=True)
labels = list(df.index)

#Plot the frequencies as barplots
values = df.iloc[:, 0]
values = np.asarray(values)
x = np.arange(len(labels))  # the label locations
import matplotlib.pyplot as plt
f, ax= plt.subplots()
ax.bar(x, values, tick_label = labels, color = "brown")
ax.set_ylabel('Number of IGRs', size = 12)



#Frequency of the type of IGRs for the conserved IGRs

my_reference_IGRS_conserved = os.path.join("reference_results", "conserved_IGRs_my_reference.fasta")
flanking_genes = []
for i in SeqIO.parse(my_reference_IGRS_conserved, "fasta"):
    flanking_genes.append(i.id.split("++")[3])

# using Counter to find frequency of elements
frequency = collections.Counter(flanking_genes)

# Convert frequency to a dictionary an then to a dataframe
frequency = dict(frequency)
df = pd.DataFrame.from_dict(frequency, orient='index')
df.rename(columns={0: "flanking_genes"}, inplace=True)
labels = list(df.index)

#Plot the frequencies as barplots
values = df.iloc[:, 0]
values = np.asarray(values)
x = np.arange(len(labels))  # the label locations
import matplotlib.pyplot as plt
f, ax= plt.subplots()
ax.bar(x, values, tick_label = labels, color = "gold")
ax.set_ylabel('Number of IGRs', size = 12)





#length distribution call get_intergenic_regions without limiting the length of the IGRS
my_input_dir = "/media/usuario/lab_bioinformatica/Carlos_Caicedo-Montoya/Genome-wide_identification_sRNAS_SClavuligerus/raw_data/Genbank_group1/"
my_reference = "GCF_005519465.1_ASM551946v1_genomic.gbff"

handle = my_input_dir + my_reference
get_intergenic_regions(handle, intergenic_length=1, max_intergenic_length=10000000, output_dir="reference_results")
my_reference_IGRS_fasta = os.path.join("reference_results", my_reference.split(".")[0] + "_IGRs.fasta")

lengths=[]
for i in SeqIO.parse(my_reference_IGRS_fasta, "fasta"):
    lengths.append(len(i))


f, ax= plt.subplots(1, 2)
ax[0].boxplot(np.asarray(lengths))
ax[0].grid()
ax[0].set_ylabel("Length of IGR (bp)")

ax[1].hist(np.asarray(lengths))
ax[1].grid()
ax[1].set_xlabel("Length of IGR (bp)")
ax[1].set_ylabel("Number of IGRs")

values = np.array([0,0.25,0.5,0.75,1])

x = np.quantile(np.asarray(lengths), values)

f.savefig ("reference_results/IGRs_boxplot.pdf", dpi=1200, format = "pdf")


#determine GC content of intergenic regions
#GC content of IGRs and conserved IGRs
from Bio.SeqUtils import GC 
from Bio.Seq import Seq

nucleotide = Seq("GACTGACTTCGA") 
GC(nucleotide) 
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

 






