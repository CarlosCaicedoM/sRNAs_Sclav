#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 11:19:43 2024

@author: usuario
"""

import pandas as pd
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import squarify    # pip install squarify (algorithm for treemap)
 
#%% read CHR data
rockhopper_CHR = pd.read_table("results_rna_seq/S-clavuligerus_RNA-Seq-data_SRA/Rockhopper3/NZ_CP027858_transcripts.txt", sep="\t") 

desired_columns = ['Transcription Start', 'Translation Start', 'Translation Stop',
                   'Transcription Stop', 'Strand', 'Name', 'Synonym', 'Product']

# Filter columns that start with 'Raw'
raw_columns = [col for col in rockhopper_CHR.columns if col.startswith('Raw')]

# Combine the lists of desired columns and columns starting with 'Raw'
columns_to_include = desired_columns + raw_columns

# Create a new DataFrame with the specified columns
rockhopper_CHR_filtered = rockhopper_CHR[columns_to_include]


rockhopper_CHR_filtered_genes = rockhopper_CHR_filtered[~rockhopper_CHR_filtered['Synonym'].str.startswith('predicted RNA')]
rockhopper_CHR_filtered_pRNAs = rockhopper_CHR_filtered[rockhopper_CHR_filtered['Synonym'].str.startswith('predicted RNA')]
rockhopper_CHR_filtered_asRNA = rockhopper_CHR_filtered_pRNAs[rockhopper_CHR_filtered_pRNAs['Product'].str.startswith('antisense:')]
rockhopper_CHR_filtered_sRNA = rockhopper_CHR_filtered_pRNAs[~rockhopper_CHR_filtered_pRNAs['Product'].str.startswith('antisense:')]

rockhopper_CHR_filtered_tRNAs = rockhopper_CHR_filtered_genes[rockhopper_CHR_filtered_genes['Product'].str.startswith('tRNA-')]
rockhopper_CHR_filtered_23S = rockhopper_CHR_filtered_genes[rockhopper_CHR_filtered_genes['Product'].str.startswith('23S ribosomal RNA')]
rockhopper_CHR_filtered_16S = rockhopper_CHR_filtered_genes[rockhopper_CHR_filtered_genes['Product'].str.startswith('16S ribosomal RNA')]
rockhopper_CHR_filtered_5S = rockhopper_CHR_filtered_genes[rockhopper_CHR_filtered_genes['Product'].str.startswith('5S ribosomal RNA')]

rockhopper_CHR_filtered_rRNA = pd.concat([rockhopper_CHR_filtered_23S, 
                                       rockhopper_CHR_filtered_16S,
                                       rockhopper_CHR_filtered_5S], ignore_index=True)

del rockhopper_CHR

#%% plasmid PLASMID 

rockhopper_plasmid = pd.read_table("results_rna_seq/S-clavuligerus_RNA-Seq-data_SRA/Rockhopper3/NZ_CP027859_transcripts.txt", sep="\t") 


desired_columns = ['Transcription Start', 'Translation Start', 'Translation Stop',
                   'Transcription Stop', 'Strand', 'Name', 'Synonym', 'Product']

# Filter columns that start with 'Raw'
raw_columns = [col for col in rockhopper_plasmid.columns if col.startswith('Raw')]

# Combine the lists of desired columns and columns starting with 'Raw'
columns_to_include = desired_columns + raw_columns

# Create a new DataFrame with the specified columns
rockhopper_plasmid_filtered = rockhopper_plasmid[columns_to_include]


rockhopper_plasmid_filtered_genes = rockhopper_plasmid_filtered[~rockhopper_plasmid_filtered['Synonym'].str.startswith('predicted RNA')]
rockhopper_plasmid_filtered_pRNAs = rockhopper_plasmid_filtered[rockhopper_plasmid_filtered['Synonym'].str.startswith('predicted RNA')]
rockhopper_plasmid_filtered_asRNA = rockhopper_plasmid_filtered_pRNAs[rockhopper_plasmid_filtered_pRNAs['Product'].str.startswith('antisense:')]
rockhopper_plasmid_filtered_sRNA = rockhopper_plasmid_filtered_pRNAs[~rockhopper_plasmid_filtered_pRNAs['Product'].str.startswith('antisense:')]

rockhopper_plasmid_filtered_tRNAs = rockhopper_plasmid_filtered_genes[rockhopper_plasmid_filtered_genes['Product'].str.startswith('tRNA-')]
rockhopper_plasmid_filtered_23S = rockhopper_plasmid_filtered_genes[rockhopper_plasmid_filtered_genes['Product'].str.startswith('23S ribosomal RNA')]
rockhopper_plasmid_filtered_16S = rockhopper_plasmid_filtered_genes[rockhopper_plasmid_filtered_genes['Product'].str.startswith('16S ribosomal RNA')]
rockhopper_plasmid_filtered_5S = rockhopper_plasmid_filtered_genes[rockhopper_plasmid_filtered_genes['Product'].str.startswith('5S ribosomal RNA')]

rockhopper_plasmid_filtered_rRNA = pd.concat([rockhopper_plasmid_filtered_23S, 
                                       rockhopper_plasmid_filtered_16S,
                                       rockhopper_plasmid_filtered_5S], ignore_index=True)

del rockhopper_plasmid

#%% rockhopper all

rockhopper_all_filtered_genes = pd.concat([rockhopper_CHR_filtered_genes, 
                                          rockhopper_plasmid_filtered_genes])

rockhopper_all_filtered_genes.set_index('Synonym', inplace=True)

#open GFF file
file_gff_all = os.path.join(os.getcwd(), "raw_data", 
                            "S-clavuligerus_ATCC27064_annotations", 
                            "GCF_005519465.1_ASM551946v1_genomic.gff")

gff_all = pd.read_csv(file_gff_all, comment='#', 
                            sep = "\t", header=None)

gff_columns = ["seqname", "source", "feature", "start", "end", "score", 
                   "strand", "frame", "attribute" ]
gff_all.columns = gff_columns

gff_all_refseq = gff_all[gff_all['source'].str.startswith('RefSeq')]

#open tsv with annotations
file_tsv_all = os.path.join(os.getcwd(), "raw_data", 
                            "S-clavuligerus_ATCC27064_annotations", 
                            "GCF_005519465.1_ASM551946v1_genomic.tsv")

tsv_annotations = pd.read_table(file_tsv_all, sep = "\t")
tsv_annotations.set_index('Locus tag', inplace=True)

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
rockhopper_all_genes_plus['complete transcript'] = 'non-complete'  # Initialize the column with "non-complete"
for index in rockhopper_all_genes_plus.index:
    if rockhopper_all_genes_plus.loc[index, 'Start'] <= tsv_annotations_plus.loc[index, 'Begin'] and rockhopper_all_genes_plus.loc[index, 'Stop'] >= tsv_annotations_plus.loc[index, 'End']:
        rockhopper_all_genes_plus.loc[index, 'complete transcript'] = 'complete'  # Update to "complete" if the condition is met

#determine if the assembly cover the annotated genes minus strand
rockhopper_all_genes_minus_formated['complete transcript'] = 'non-complete'  # Initialize the column with "non-complete"
for index in rockhopper_all_genes_minus_formated.index:
    if rockhopper_all_genes_minus_formated.loc[index, 'Start'] <= tsv_annotations_minus.loc[index, 'Begin'] and rockhopper_all_genes_minus_formated.loc[index, 'Stop'] >= tsv_annotations_minus.loc[index, 'End']:
        rockhopper_all_genes_minus_formated.loc[index, 'complete transcript'] = 'complete'  # Update to "complete" if the condition is met


#%% analyze aligned reads

aligned_reads_lines = []

# path to the file
file_path = 'results_rna_seq/S-clavuligerus_RNA-Seq-data_SRA/Rockhopper3/summary.txt' 
with open(file_path, 'r') as file:
    for line in file:
        # Verificar si la línea comienza con "Successfully aligned reads:"
        if line.startswith('Successfully aligned reads:'):
            # Agregar la línea a la lista
            aligned_reads_lines.append(line.strip())  # Strip para eliminar espacios en blanco al principio y al final

aligned_reads_lines_chr = []
aligned_reads_lines_plasmid = []
for line in aligned_reads_lines:
    if line.endswith("chromosome)"):
        aligned_reads_lines_chr.append(line)
    elif line.endswith("pCLA1)"):
        aligned_reads_lines_plasmid.append(line)

def extract_percentage(line):
    # Use a regular expression to find the number before the '%' symbol
    match = re.search(r'(\d+)%', line)
    if match:
        return int(match.group(1))  # Convert the value to an integer
    else:
        return None  # Return None if the pattern is not found

# Apply the function to each element in the lists
percentage_chr = [extract_percentage(line) for line in aligned_reads_lines_chr]
percentage_plasmid = [extract_percentage(line) for line in aligned_reads_lines_plasmid]

percentage_chr = np.asarray(percentage_chr)
percentage_plasmid = np.asarray(percentage_plasmid)

percentage_total = percentage_chr+percentage_plasmid

average_aligned_reads = np.mean(percentage_total)


#%% Analyze expression

#all genes  expression
rockhopper_all_filtered_genes['Number of samples in which expressed'] = (rockhopper_all_filtered_genes.filter(like='Raw').ne(0).sum(axis=1))
rockhopper_all_filtered_genes = rockhopper_all_filtered_genes.drop(['Name'], axis=1)
cumulative_distribution_all_genes  = rockhopper_all_filtered_genes['Number of samples in which expressed'].value_counts(normalize=True).sort_index().cumsum()

## CDS expresion
tsv_CDS = tsv_annotations[tsv_annotations['Gene Type'].str.startswith('protein-coding')]
rockhopper_all_filtered_CDS = rockhopper_all_filtered_genes.join(tsv_CDS, how='inner')
rockhopper_all_filtered_CDS = rockhopper_all_filtered_CDS.drop(tsv_CDS.columns, axis=1)
cumulative_distribution_CDS = rockhopper_all_filtered_CDS['Number of samples in which expressed'].value_counts(normalize=True).sort_index().cumsum()

## tRNAs expression
rockhopper_all_filtered_tRNAs = pd.concat([rockhopper_CHR_filtered_tRNAs, 
                                           rockhopper_plasmid_filtered_tRNAs])
rockhopper_all_filtered_tRNAs['Number of samples in which expressed'] = (rockhopper_all_filtered_tRNAs.filter(like='Raw').ne(0).sum(axis=1))
cumulative_distribution_all_tRNAs  = rockhopper_all_filtered_tRNAs['Number of samples in which expressed'].value_counts(normalize=True).sort_index().cumsum()


## rRNAs expression
rockhopper_all_filtered_rRNAs = pd.concat([rockhopper_CHR_filtered_rRNA, 
                                           rockhopper_plasmid_filtered_rRNA])
rockhopper_all_filtered_rRNAs['Number of samples in which expressed'] = (rockhopper_all_filtered_rRNAs.filter(like='Raw').ne(0).sum(axis=1))
cumulative_distribution_all_rRNAs  = rockhopper_all_filtered_rRNAs['Number of samples in which expressed'].value_counts(normalize=True).sort_index().cumsum()

## as expression
rockhopper_all_filtered_asRNAs = pd.concat([rockhopper_CHR_filtered_asRNA, 
                                           rockhopper_plasmid_filtered_asRNA])
rockhopper_all_filtered_asRNAs['Number of samples in which expressed'] = (rockhopper_all_filtered_asRNAs.filter(like='Raw').ne(0).sum(axis=1))
cumulative_distribution_all_asRNAs  = rockhopper_all_filtered_asRNAs['Number of samples in which expressed'].value_counts(normalize=True).sort_index().cumsum()

#novel transcripts expression (sRNAs)
rockhopper_all_filtered_sRNAs = pd.concat([rockhopper_CHR_filtered_sRNA, 
                                           rockhopper_plasmid_filtered_sRNA])
rockhopper_all_filtered_sRNAs['Number of samples in which expressed'] = (rockhopper_all_filtered_sRNAs.filter(like='Raw').ne(0).sum(axis=1))
cumulative_distribution_all_sRNAs  = rockhopper_all_filtered_asRNAs['Number of samples in which expressed'].value_counts(normalize=True).sort_index().cumsum()




# Plot the cumulative distribution
f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)

ax1.plot(cumulative_distribution_all_genes.index, 
         cumulative_distribution_all_genes.values, 
         marker = "+", label = "All genes")
ax1.set_xlabel('Number of samples')
ax1.set_ylabel('Cumulative Distribution')
ax1.grid()
ax1.legend()

ax2.plot(cumulative_distribution_CDS.index, 
         cumulative_distribution_CDS.values, 
         marker='+', label  = "CDS")
ax2.set_xlabel('Number of samples')
ax2.grid()
ax2.legend()

ax3.plot(cumulative_distribution_all_tRNAs.index, 
         cumulative_distribution_all_tRNAs.values, 
         marker='+', label = "tRNA")
ax3.set_xlabel('Number of samples')
ax3.grid()
ax3.legend()

ax4.plot(cumulative_distribution_all_rRNAs.index, 
         cumulative_distribution_all_rRNAs.values, 
         marker='+', label = "rRNA")
ax4.set_xlabel('Number of samples')
ax4.set_ylabel('Cumulative Distribution')
ax4.grid()
ax4.legend()

ax5.plot(cumulative_distribution_all_asRNAs.index, 
         cumulative_distribution_all_asRNAs.values, 
         marker='+', label = "antisense")
ax5.set_xlabel('Number of samples')
ax5.grid()
ax5.legend()

ax6.plot(cumulative_distribution_all_sRNAs.index, 
         cumulative_distribution_all_sRNAs.values,
         marker='+', label = "Novel transcripts")
ax6.set_xlabel('Number of samples')
ax6.grid()
ax6.legend()

#%% TPM calculation

#reads per kilobase (RPK)

rockhopper_all_genes_plus_minus = pd.concat([rockhopper_all_genes_plus, 
                                           rockhopper_all_genes_minus_formated])

def TPM(df):
    # Calculate the difference in kilobases
    df['Diff'] = (df['Stop'] - df['Start']) / 1000
    
    # Calculate RPK values
    for col in df.columns:
        if col.startswith('Raw Counts'):
            new_col_name = col.replace('Raw Counts', 'RPK')
            df[new_col_name] = df[col] / df['Diff']

    # Drop the 'Diff' column
    df.drop('Diff', axis=1, inplace=True)

    # Function to calculate the scaling factor
    def calculate_scaling_factor(row):
        rpm_cols = [col for col in row.index if col.startswith('RPK')]
        total_rpm = sum(row[rpm_cols])
        scaling_factor = total_rpm / 1000000
        return scaling_factor

    # Calculate the scaling factor
    df['scaling_factor'] = df.apply(calculate_scaling_factor, axis=1)

    # Calculate TPM values
    for col in df.columns:
        if col.startswith('RPK'):
            new_col_name = col.replace('RPK', 'TPM')
            df[new_col_name] = df[col] / df['scaling_factor']

    # Calculate the maximum TPM for each row
    df['maximum_TPM'] = df.apply(lambda row: row[[col for col in df.columns if col.startswith('TPM')]].max(), axis=1)

    return df


def define_transcript_start_stop(df):
    # Define the transcript start 
    df['Start'] = df.apply(lambda row: row['Transcription Start'] if pd.notna(row['Transcription Start']) else row['Translation Start'], axis=1)
    columns = df.columns.tolist()
    start_column = columns.pop(columns.index('Start'))
    columns.insert(6, start_column)
    df = df.reindex(columns=columns)

    # Define the transcript stop 
    df['Stop'] = df.apply(lambda row: row['Transcription Stop'] if pd.notna(row['Transcription Stop']) else row['Translation Stop'], axis=1)
    columns = df.columns.tolist()
    stop_column = columns.pop(columns.index('Stop'))
    columns.insert(7, stop_column)
    df = df.reindex(columns=columns)

    return df

#define transcript start stop
rockhopper_all_filtered_CDS = define_transcript_start_stop(rockhopper_all_filtered_CDS)
rockhopper_all_filtered_tRNAs = define_transcript_start_stop(rockhopper_all_filtered_tRNAs)
rockhopper_all_filtered_rRNAs = define_transcript_start_stop(rockhopper_all_filtered_rRNAs)
rockhopper_all_filtered_asRNAs = define_transcript_start_stop(rockhopper_all_filtered_asRNAs)
rockhopper_all_filtered_sRNAs = define_transcript_start_stop(rockhopper_all_filtered_sRNAs)




#calculate TPMs
#all genes
rockhopper_all_genes_plus_minus = TPM(rockhopper_all_genes_plus_minus)

#CDS
rockhopper_all_filtered_CDS = TPM(rockhopper_all_filtered_CDS)

#tRNAs
rockhopper_all_filtered_tRNAs = TPM(rockhopper_all_filtered_tRNAs)

#rRNAs
rockhopper_all_filtered_rRNAs = TPM(rockhopper_all_filtered_rRNAs)

#asRNAs
rockhopper_all_filtered_asRNAs = TPM(rockhopper_all_filtered_asRNAs)

#sRNAs
rockhopper_all_filtered_sRNAs = TPM(rockhopper_all_filtered_sRNAs)


# Dibujar el histograma de la columna "maximum_TPM"
f, ((ax1, ax2, ax3), (ax4, ax5, ax6) ) = plt.subplots(2, 3)
for ax in (ax1, ax2, ax3, ax4, ax5, ax6):
    ax.tick_params(axis='both', labelsize=8)  # Cambiar el tamaño de la fuente a 8 (o el tamaño deseado)

ax1.hist(rockhopper_all_genes_plus_minus['maximum_TPM'],
         bins=100, edgecolor='black', label = "All genes")
ax1.set_xlabel('Maximum TPM')
ax1.set_ylabel('Transcripts Frequency')    
ax1.grid()
ax1.legend()

ax2.hist(rockhopper_all_filtered_CDS['maximum_TPM'],
         bins=100, edgecolor='black', label = "CDS", color = "orange")
ax2.set_xlabel('Maximum TPM')
ax2.grid()
ax2.legend()

ax3.hist(rockhopper_all_filtered_tRNAs['maximum_TPM'],
         bins=100, edgecolor='black', label = "tRNA", color = "darkcyan")
ax3.set_xlabel('Maximum TPM')
ax3.grid()
ax3.legend()

ax4.hist(rockhopper_all_filtered_rRNAs['maximum_TPM'],
         bins=100, edgecolor='black', label = "rRNA", color = "lawngreen")
ax4.set_xlabel('Maximum TPM')
ax4.set_ylabel('Transcripts Frequency')    
ax4.grid()
ax4.legend()

ax5.hist(rockhopper_all_filtered_asRNAs['maximum_TPM'],
         bins=100, edgecolor='black', label = "Antisense", color = "darkred")
ax5.set_xlabel('Maximum TPM')
ax5.grid()
ax5.legend()

ax6.hist(rockhopper_all_filtered_sRNAs['maximum_TPM'],
         bins=100, edgecolor='black', label = "Novel transcripts", color = "indigo")
ax6.set_xlabel('Maximum TPM')
ax6.grid()
ax6.legend()


#%% operons prediction


def read_results_operons(path):
    operons_predictions = pd.read_table(path,  sep="\t", header=0)
    return operons_predictions

chromosome = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/results_rna_seq/S-clavuligerus_RNA-Seq-data_SRA/Rockhopper3/NZ_CP027858_operons.txt" 
plasmid = "/home/usuario/Documentos/Carlos_PhD/sRNAS_Sclav/results_rna_seq/S-clavuligerus_RNA-Seq-data_SRA/Rockhopper3/NZ_CP027859_operons.txt"

operons_CHR = read_results_operons(chromosome)
operons_plasmid = read_results_operons(plasmid)

total_genes_operons_plasmid = operons_plasmid['Number of Genes'].sum()
total_genes_operons_CHR = operons_CHR['Number of Genes'].sum()


gene_pairs_operons_CHR = operons_CHR[operons_CHR['Number of Genes'] == 2]
gene_pairs_operons_plasmid = operons_plasmid[operons_plasmid['Number of Genes'] == 2]
multi_gene_operons_CHR = operons_CHR[operons_CHR['Number of Genes'] != 2]
multi_gene_operons_plasmid = operons_plasmid[operons_plasmid['Number of Genes'] != 2]



total_gene_pairs_operons_CHR = gene_pairs_operons_CHR['Number of Genes'].sum()
total_gene_pairs_operons_plasmid = gene_pairs_operons_plasmid['Number of Genes'].sum()
 
total_multi_gene_operons_CHR = multi_gene_operons_CHR['Number of Genes'].sum()
total_multi_gene_operons_plasmid = multi_gene_operons_plasmid['Number of Genes'].sum()




# Create a data frame with fake data
total_operons = pd.DataFrame({'totals':[total_gene_pairs_operons_CHR,
                                   total_multi_gene_operons_CHR,
                                   total_gene_pairs_operons_plasmid,
                                   total_multi_gene_operons_plasmid], 
                         'group':["gene-pairs chromosome", 
                                  "multi-gene operons chromosome", 
                                  "gene-pairs plasmid",
                                  "multi-gene operons plasmid"] })


# plot it
squarify.plot(sizes=total_operons['totals'], label=total_operons['group'])
plt.axis('off')
plt.show()

# Dibujar el histograma de la columna "maximum_TPM"
f, (ax1, ax2) = plt.subplots(1, 2)

ax1.hist(operons_CHR['Number of Genes'],
         bins=20, edgecolor='black', label = "Chromosome")
ax1.set_xlabel('Operon size')
ax1.set_ylabel('Frequency')    
ax1.grid()
ax1.legend()

ax2.hist(operons_plasmid['Number of Genes'],
         bins=20, edgecolor='black', 
         label = "Plasmid", color = "darkorange")
ax2.set_xlabel('Operon size')
ax2.grid()
ax2.legend()