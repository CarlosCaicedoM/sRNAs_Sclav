#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 11:19:43 2024

@author: usuario
"""

import pandas as pd
import os


#CHR
rockhopper_results_CHR = pd.read_table("results_rna_seq/S-clavuligerus_RNA-Seq-data_SRA/Rockhopper3/NZ_CP027858_transcripts.txt", sep="\t") 

desired_columns = ['Transcription Start', 'Translation Start', 'Translation Stop',
                   'Transcription Stop', 'Strand', 'Name', 'Synonym', 'Product']

# Filter columns that start with 'Raw'
raw_columns = [col for col in rockhopper_results_CHR.columns if col.startswith('Raw')]

# Combine the lists of desired columns and columns starting with 'Raw'
columns_to_include = desired_columns + raw_columns

# Create a new DataFrame with the specified columns
rockhopper_results_CHR_filtered = rockhopper_results_CHR[columns_to_include]


rockhopper_results_CHR_filtered_genes = rockhopper_results_CHR_filtered[~rockhopper_results_CHR_filtered['Synonym'].str.startswith('predicted RNA')]
rockhopper_results_CHR_filtered_pRNAs = rockhopper_results_CHR_filtered[rockhopper_results_CHR_filtered['Synonym'].str.startswith('predicted RNA')]
rockhopper_results_CHR_filtered_asRNA = rockhopper_results_CHR_filtered_pRNAs[rockhopper_results_CHR_filtered_pRNAs['Product'].str.startswith('antisense:')]
rockhopper_results_CHR_filtered_sRNA = rockhopper_results_CHR_filtered_pRNAs[~rockhopper_results_CHR_filtered_pRNAs['Product'].str.startswith('antisense:')]

del rockhopper_results_CHR

#PLASMID 

rockhopper_results_plasmid = pd.read_table("results_rna_seq/S-clavuligerus_RNA-Seq-data_SRA/Rockhopper3/NZ_CP027859_transcripts.txt", sep="\t") 

# Filter columns that start with 'Raw'
raw_columns = [col for col in rockhopper_results_plasmid.columns if col.startswith('Raw')]

# Combine the lists of desired columns and columns starting with 'Raw'
columns_to_include = desired_columns + raw_columns

# Create a new DataFrame with the specified columns
rockhopper_results_plasmid_filtered = rockhopper_results_plasmid[columns_to_include]

rockhopper_results_plasmid_filtered_genes = rockhopper_results_plasmid_filtered[~rockhopper_results_plasmid_filtered['Synonym'].str.startswith('predicted RNA')]
rockhopper_results_plasmid_filtered_pRNAs = rockhopper_results_plasmid_filtered[rockhopper_results_plasmid_filtered['Synonym'].str.startswith('predicted RNA')]
rockhopper_results_plasmid_filtered_asRNA = rockhopper_results_plasmid_filtered_pRNAs[rockhopper_results_plasmid_filtered_pRNAs['Product'].str.startswith('antisense:')]
rockhopper_results_plasmid_filtered_sRNA = rockhopper_results_plasmid_filtered_pRNAs[~rockhopper_results_plasmid_filtered_pRNAs['Product'].str.startswith('antisense:')]



del rockhopper_results_plasmid

#rockhopper all

rockhopper_all_filtered_genes = pd.concat([rockhopper_results_CHR_filtered_genes, 
                                          rockhopper_results_plasmid_filtered_genes])

#open GFF file
file_gff_all = os.path.join(os.getcwd(), "raw_data", 
                            "S-clavuligerus_ATCC27064_annotations", 
                            "GCF_005519465.1_ASM551946v1_genomic.gff")

gff_all = pd.read_csv(file_gff_all, comment='#', 
                            sep = "\t", header=None)

gff_columns = ["seqname", "source", "feature", "start", "end", "score", 
                   "strand", "frame", "attribute" ]
gff_all.columns = gff_columns

gff_all_refseq = gff_all[gff_all['source'].str.startswith('Refseq')]

#open tsv with annotations
file_tsv_all = os.path.join(os.getcwd(), "raw_data", 
                            "S-clavuligerus_ATCC27064_annotations", 
                            "GCF_005519465.1_ASM551946v1_genomic.tsv")

tsv_annotations = pd.read_table(file_tsv_all, sep = "\t")

tsv_annotations_plus = tsv_annotations[tsv_annotations['Orientation'].str.startswith('plus')]
tsv_annotations_minus = tsv_annotations[tsv_annotations['Orientation'].str.startswith('minus')]


rockhopper_all_genes_plus = rockhopper_all_filtered_genes[rockhopper_all_filtered_genes['Strand'].str.startswith('+')]
rockhopper_all_genes_minus = rockhopper_all_filtered_genes[rockhopper_all_filtered_genes['Strand'].str.startswith('-')]


#Transform minus annotations to be similar to the ones in tsv file
rockhopper_all_genes_minus_formated = rockhopper_all_genes_minus.copy() 
column_name_mapping = {col: col.replace("Start", "Stop") if col.endswith("Start") else col.replace("Stop", "Start") for col in rockhopper_all_genes_minus_formated.columns}
rockhopper_all_genes_minus_formated = rockhopper_all_genes_minus_formated.rename(columns=column_name_mapping)

#define the transcript start for the plus strand
rockhopper_all_genes_plus['Start'] = rockhopper_all_genes_plus.apply(lambda row: row['Transcription Start'] if pd.notna(row['Transcription Start']) else row['Translation Start'], axis=1)
columns = rockhopper_all_genes_plus.columns.tolist()
start_column = columns.pop(columns.index('Start'))
columns.insert(6, start_column)
rockhopper_all_genes_plus = rockhopper_all_genes_plus.reindex(columns=columns)

#define the transcript stop for the plus strand
rockhopper_all_genes_plus['Stop'] = rockhopper_all_genes_plus.apply(lambda row: row['Transcription Stop'] if pd.notna(row['Transcription Stop']) else row['Translation Stop'], axis=1)
columns = rockhopper_all_genes_plus.columns.tolist()
start_column = columns.pop(columns.index('Stop'))
columns.insert(7, start_column)
rockhopper_all_genes_plus = rockhopper_all_genes_plus.reindex(columns=columns)


#define the transcript start for the minus strand
rockhopper_all_genes_minus_formated['Start'] = rockhopper_all_genes_minus_formated.apply(lambda row: row['Transcription Start'] if pd.notna(row['Transcription Start']) else row['Translation Start'], axis=1)
columns = rockhopper_all_genes_minus_formated.columns.tolist()
start_column = columns.pop(columns.index('Start'))
columns.insert(6, start_column)
rockhopper_all_genes_minus_formated = rockhopper_all_genes_minus_formated.reindex(columns=columns)

#define the transcript stop for the minus strand
rockhopper_all_genes_minus_formated['Stop'] = rockhopper_all_genes_minus_formated.apply(lambda row: row['Transcription Stop'] if pd.notna(row['Transcription Stop']) else row['Translation Stop'], axis=1)
columns = rockhopper_all_genes_minus_formated.columns.tolist()
start_column = columns.pop(columns.index('Stop'))
columns.insert(7, start_column)
rockhopper_all_genes_minus_formated = rockhopper_all_genes_minus_formated.reindex(columns=columns)

#determine if the assembly cover the annotated genes plus strand
rockhopper_all_genes_plus['complete transcript'] = rockhopper_all_genes_plus.apply(lambda row: 'complete' if (row['Start'] <= tsv_annotations_plus['Begin']).all() and (row['Stop'] >= tsv_annotations_plus['End']).all() else 'non-complete', axis=1)







