#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 23:26:58 2024

@author: usuario
"""

#%% Read and format the data
import pandas as pd
import os
import matplotlib.pyplot as plt


#Add data of promoters and terminators detection
all_IGRs_term_prom = pd.read_table("results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/reference_results/all_IGRs_term_prom.csv", sep=",")
all_IGRs_term_prom = all_IGRs_term_prom.rename(columns={all_IGRs_term_prom.columns[0]: 'header'})
all_IGRs_term_prom['Promoter1'] = all_IGRs_term_prom.apply(lambda row: (row['start'], row['end'], row['strand'], row['sequence']), axis=1)
all_IGRs_term_prom.drop(columns=['start', 'end', 'strand', 'sequence'], inplace=True)
columns = [all_IGRs_term_prom.columns[0]] + ['Promoter1'] + [col for col in all_IGRs_term_prom if col != 'Promoter1' and col != all_IGRs_term_prom.columns[0]]
all_IGRs_term_prom = all_IGRs_term_prom[columns]
all_IGRs_term_prom.set_index('header', inplace=True)

#add data about conservation and RNAz predictions

conserved_IGRs = pd.read_table("results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/reference_results/putative_sRNAs_rfam_blast_infernal_CP_NOG.csv", sep=",")
conserved_IGRs = conserved_IGRs.rename(columns={conserved_IGRs.columns[0]: 'header'})

columns_to_extract = ['header', 'RNAz', 'Strand1', 'RNA_Class_probability1', 
                      'Strand2', 'RNA_Class_probability2', 'number_of_genomes']
conserved_IGRs = conserved_IGRs[columns_to_extract]
conserved_IGRs.set_index('header', inplace=True)

IGRS_prom_term_RNAz = pd.concat([all_IGRs_term_prom, 
                                                conserved_IGRs], axis = 1)


# Add data of coding potential

coding_potential = pd.read_table("results_sRNA_clade2_pangenome_gene_nontrimmed_50-600/CPC2/GCF_005519465_ALL-IGRS_cpc2.txt", sep="\t")
coding_potential_reordered = coding_potential[['#ID', 'label', 'peptide_length', 'Fickett_score', 'pI', 'ORF_integrity', 'coding_probability']]

coding_potential_reordered.set_index('#ID', inplace=True)


IGRS_prom_term_RNAz_CP = pd.concat([IGRS_prom_term_RNAz, 
                                                coding_potential_reordered], axis = 1)

#Add blast annotations

blast_all_IGRs = pd.read_table(os.path.join("results_sRNA_clade2_pangenome_gene_nontrimmed_50-600", 
                                          "reference_results", 
                                          "all_IGRs_blastn_out.txt"), sep='\t', header=None)

blast_all_IGRs.columns = ["query", "subject", "identity", "coverage",
                   "qlength", "slength", "alength",
                   "bitscore", "E-value"]
blast_all_IGRs.set_index("query", inplace=True)

IGRS_prom_term_RNAz_CP_B = pd.concat([IGRS_prom_term_RNAz_CP, 
                                                blast_all_IGRs], axis = 1)


#Add Infernal annotations

rfam_results = os.path.join("results_sRNA_clade2_pangenome_gene_nontrimmed_50-600",
                                          "Infernal", 
                                          "Sclav-all_IGRs.tblout")

infernal_annotations = pd.read_csv(rfam_results, comment='#', header=None)
infernal_annotations = infernal_annotations[0].str.split(expand=True)
infernal_annotations['description_of_target'] = infernal_annotations.apply(lambda row: ' '.join(str(val) for val in row[17:-1] if val is not None), axis=1)
infernal_annotations_filtered = infernal_annotations.iloc[:, [-1, 1,2,7,8,9,14,15]]

infernal_columns = ["target_name",  "accession", "query_name", "seq_from",
                    "seq_to", "strand", "score", "E-value"]
infernal_annotations_filtered.columns = infernal_columns
infernal_annotations_filtered.index.duplicated()
infernal_annotations_filtered = infernal_annotations_filtered.drop_duplicates(subset=['query_name'])
infernal_annotations_filtered = infernal_annotations_filtered.set_index('query_name')

IGRS_prom_term_RNAz_CP_B_INF = pd.concat([IGRS_prom_term_RNAz_CP_B, 
                                                infernal_annotations_filtered], axis = 1)



#Add Rockhopper results

##CHR
rockhopper_CHR = pd.read_table("results_rna_seq/S-clavuligerus_RNA-Seq-data_SRA/Rockhopper3/NZ_CP027858_transcripts.txt", sep="\t") 
desired_columns = ['Transcription Start', 'Translation Start', 'Translation Stop',
                   'Transcription Stop', 'Strand', 'Name', 'Synonym', 'Product']
# Filter columns that start with 'Raw'
raw_columns = [col for col in rockhopper_CHR.columns if col.startswith('Raw')]

# Combine the lists of desired columns and columns starting with 'Raw'
columns_to_include = desired_columns + raw_columns
# Create a new DataFrame with the specified columns
rockhopper_CHR_filtered = rockhopper_CHR[columns_to_include]
rockhopper_CHR_filtered_nonAS = rockhopper_CHR_filtered[~rockhopper_CHR_filtered['Product'].str.startswith('antisense:')]


##PLASMID
rockhopper_plasmid = pd.read_table("results_rna_seq/S-clavuligerus_RNA-Seq-data_SRA/Rockhopper3/NZ_CP027859_transcripts.txt", sep="\t") 
desired_columns = ['Transcription Start', 'Translation Start', 'Translation Stop',
                   'Transcription Stop', 'Strand', 'Name', 'Synonym', 'Product']
# Filter columns that start with 'Raw'
raw_columns = [col for col in rockhopper_CHR.columns if col.startswith('Raw')]

# Combine the lists of desired columns and columns starting with 'Raw'
columns_to_include = desired_columns + raw_columns
# Create a new DataFrame with the specified columns
rockhopper_plasmid_filtered = rockhopper_plasmid[columns_to_include]
rockhopper_plasmid_filtered_nonAS = rockhopper_plasmid_filtered[~rockhopper_plasmid_filtered['Product'].str.startswith('antisense:')]



del rockhopper_CHR
del rockhopper_plasmid

rockhopper_all_filtered_nonAS = pd.concat([rockhopper_CHR_filtered_nonAS, 
                                          rockhopper_plasmid_filtered_nonAS])

#%% 


df = rockhopper_all_filtered_nonAS
# Reset the index of the DataFrame
df = df.reset_index(drop=True)

# Create a new column in the DataFrame to store the flanking entries for "predicted RNA"
df['Flanking'] = ''

# Iterate over the rows of the DataFrame
for index, row in df.iterrows():
    # Check if the row contains "predicted RNA" in the "Synonym" column
    if row['Synonym'] == 'predicted RNA':
        # Find the previous entry that is not "predicted RNA"
        previous_entry = ''
        for i in range(index - 1, -1, -1):
            if df.at[i, 'Synonym'] != 'predicted RNA':
                previous_entry = df.at[i, 'Synonym']
                break
        
        # Find the next entry that is not "predicted RNA"
        next_entry = ''
        for i in range(index + 1, len(df)):
            if df.at[i, 'Synonym'] != 'predicted RNA':
                next_entry = df.at[i, 'Synonym']
                break
        
        # Store the flanking entries for "predicted RNA" in the new "Flanking" column
        df.at[index, 'Flanking'] = f'{previous_entry}-{next_entry}'

# Move the "Flanking" column to position 8
flanking_column = df.pop('Flanking')
df.insert(8, 'Flanking', flanking_column)



IGR_sRNAS_rockhopper = df[df['Product'] == '-']

columns_to_drop = ['Translation Start', 'Translation Stop']

# Create a new DataFrame without the specified columns
IGR_sRNAS_rockhopper = IGR_sRNAS_rockhopper.drop(columns=columns_to_drop)

columns = list(IGR_sRNAS_rockhopper.columns)
modified_columns = [string.replace(" ", "_") for string in columns]

IGR_sRNAS_rockhopper.columns = modified_columns


# Get a list of columns starting with 'Raw'
raw_columns2 = [col for col in IGR_sRNAS_rockhopper.columns if col.startswith('Raw')]

# Filter rows where the sum of values in 'Raw' columns is greater than 10
filtered_sRNAs = IGR_sRNAS_rockhopper[IGR_sRNAS_rockhopper[raw_columns2].gt(0).sum(axis=1) > 10]
final_sRNAs = filtered_sRNAs[(filtered_sRNAs[raw_columns2] > 20).sum(axis=1) >= 20]
# Counter for numbering "Predicted RNA" entries
counter = 1

# Iterate over the rows of the DataFrame
for index, row in final_sRNAs.iterrows():
    # Check if the entry in the "Synonym" column is "Predicted RNA"
    if row['Synonym'] == 'predicted RNA':
        # Generate the new value for the "Synonym" entry
        new_value = f'predicted_RNA_{counter}'
        # Update the entry in the DataFrame
        final_sRNAs.at[index, 'Synonym'] = new_value
        # Increment the counter for the next entry
        counter += 1    

#%% Compare bioinformatics predictions and Rockhopper
bioinformatics_predictions = IGRS_prom_term_RNAz_CP_B_INF
bioinformatics_predictions.reset_index(inplace=True)
bioinformatics_predictions.rename(columns={'index': 'header'}, inplace=True)

headers_IGRs_conserved = []
flanking_genes = []
coordinates = []
orientation  = []
for header in bioinformatics_predictions['header']:
    headers_IGRs_conserved.append(header)
    flanking_genes.append(header.split("++")[1])
    coordinates.append(header.split("++")[2])
    orientation.append(header.split("++")[3])
    
IGRs_conserved_dict2 = dict(zip(headers_IGRs_conserved, coordinates))
    
    
#Locate sRNAs in Bioinformatics predictions in the Rockhopper results
empty_lists =  [[] for _ in range(len(coordinates))]
predictions_in_Rockhopper = {headers_IGRs_conserved[i]:empty_lists[i] for i in range(len(headers_IGRs_conserved))}
    
for key, value in IGRs_conserved_dict2.items():
    flanking_genes = key.split("++")[1]
    for j in final_sRNAs.itertuples():
        if (flanking_genes == j.Flanking):
            conserved_region = (j.Transcription_Start, j.Transcription_Stop, j.Strand, j.Synonym)
            predictions_in_Rockhopper[key].append(conserved_region)

 
df2_predictions_in_Rockhopper= pd.DataFrame.from_dict(predictions_in_Rockhopper, orient="index")
new_column_labels = [f"Rockhopper_{i+1}" for i in range(len(df2_predictions_in_Rockhopper.columns))]
df2_predictions_in_Rockhopper.columns = new_column_labels

counter = 0

# Iterate over the dictionary and count elements that are not empty lists
for value_list in predictions_in_Rockhopper.values():
    if value_list:  # Check if the list has elements (i.e., not empty)
        counter += 1

verified_sRNAs  =  counter

print(verified_sRNAs)

#Concatenate bioinformatics_predictions with Rockhopper Predictions
bioinformatics_predictions.set_index('header', inplace=True)
IGRS_prom_term_RNAz_CP_B_INF_R = pd.concat([bioinformatics_predictions, 
                                                df2_predictions_in_Rockhopper], axis = 1)


#%% abrir predicciones apero GFF y comparar cuales de las predicciones
my_experiment = "S-clav_RNA-seq_sRNas_OMEGA"

file_gff_APERO_PLASMID = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "APERO", "PLASMID.gff")
file_gff_APERO_CHR = os.path.join(os.getcwd(), "results_rna_seq", my_experiment, "APERO", "CHR.gff")
file_gff_all = os.path.join(os.getcwd(), "raw_data", "S-clavuligerus_ATCC27064_annotations", "GCF_005519465.1_ASM551946v1_genomic.gff")


gff_APERO_PLASMID = pd.read_table(file_gff_APERO_PLASMID, comment='#', 
                                  sep = "\t", header=None)
gff_APERO_CHR = pd.read_csv(file_gff_APERO_CHR, comment='#', 
                            sep = "\t", header=None)

gff_all = pd.read_csv(file_gff_all, comment='#', 
                            sep = "\t", header=None)
gff_columns = ["seqname", "source", "feature", "start", "end", "score", 
                   "strand", "frame", "attribute" ]
gff_all.columns = gff_columns


gff_all_refseq = gff_all[gff_all['source'] == 'RefSeq']

gff_all_refseq_coding = gff_all_refseq[~gff_all_refseq['feature'].str.contains('region|direct_repeat')]

gff_all_refseq_coding['name'] = gff_all_refseq_coding['attribute'].str.extract(r'ID=gene-([^;]+)')


gff_APERO_PLASMID.columns = gff_columns
gff_APERO_CHR.columns = gff_columns


APERO_Sclav = pd.concat([gff_APERO_CHR, gff_APERO_PLASMID])
APERO_Sclav = APERO_Sclav.reset_index(drop=True)

# Merge the DataFrames based on common columns 'source', 'start', and 'end'
APERO_Sclav2 = pd.merge(APERO_Sclav, gff_all_refseq_coding[['source', 'start', 'end', 'name']], on=['source', 'start', 'end'], how='left')
APERO_Sclav2['name'] = APERO_Sclav2['name'].fillna('putative_sRNA')
# Create a new column in the DataFrame to store the flanking entries for "predicted RNA"
APERO_Sclav2['Flanking'] = ''

# Iterate over the rows of the DataFrame
for index, row in APERO_Sclav2.iterrows():
    # Check if the row contains "predicted RNA" in the "Synonym" column
    if row['name'] == 'putative_sRNA':
        # Find the previous entry that is not "predicted RNA"
        previous_entry = ''
        for i in range(index - 1, -1, -1):
            if APERO_Sclav2.at[i, 'name'] != 'putative_sRNA':
                previous_entry = APERO_Sclav2.at[i, 'name']
                break
        
        # Find the next entry that is not "predicted RNA"
        next_entry = ''
        for i in range(index + 1, len(APERO_Sclav2)):
            if APERO_Sclav2.at[i, 'name'] != 'putative_sRNA':
                next_entry = APERO_Sclav2.at[i, 'name']
                break
        
        # Store the flanking entries for "predicted RNA" in the new "Flanking" column
        APERO_Sclav2.at[index, 'Flanking'] = f'{previous_entry}-{next_entry}'

# Move the "Flanking" column to position 8
sRNAS_APERO = APERO_Sclav2[APERO_Sclav2['name'] == 'putative_sRNA']
sRNAS_APERO = sRNAS_APERO[sRNAS_APERO['feature'] == 'putative_sRNA']

#%% Compare bioinformatics predictions and Rockhopper

#Locate sRNAs in Bioinformatics predictions in the Rockhopper results
empty_lists2 =  [[] for _ in range(len(coordinates))]
predictions_in_APERO = {headers_IGRs_conserved[i]:empty_lists2[i] for i in range(len(headers_IGRs_conserved))}
    
for key, value in IGRs_conserved_dict2.items():
    flanking_genes = key.split("++")[1]
    for j in sRNAS_APERO.itertuples():
        if (flanking_genes == j.Flanking):
            conserved_region = (j.start, j.end, j.strand, j.name)
            predictions_in_APERO[key].append(conserved_region)

 
df2_predictions_in_APERO= pd.DataFrame.from_dict(predictions_in_APERO, orient="index")
APERO_labels = [f"APERO_{i+1}" for i in range(len(df2_predictions_in_APERO.columns))]
df2_predictions_in_APERO.columns = APERO_labels

#Concatenate bioinformatics_predictions with Rockhopper Predictions
IGRS_prom_term_RNAz_CP_B_INF_R_A = pd.concat([IGRS_prom_term_RNAz_CP_B_INF_R, 
                                                df2_predictions_in_APERO], axis = 1)




#%% Summary



apero_columns = [col for col in IGRS_prom_term_RNAz_CP_B_INF_R_A.columns if col.startswith('APERO')]
R_columns = [col for col in IGRS_prom_term_RNAz_CP_B_INF_R_A.columns if col.startswith('Rockhopper')]


IGRs_summary = IGRS_prom_term_RNAz_CP_B_INF_R_A[[ 'Promoter1', 
                                                     'Terminator1',
                                                     'number_of_genomes',
                                                     'RNAz',
                                                     'Strand1',
                                                     'RNA_Class_probability1',
                                                     'Strand2',
                                                     'RNA_Class_probability2',
                                                     'label',
                                                     'subject',
                                                     'target_name',
                                                    ]+R_columns+
                                                apero_columns].copy()


#%% See the predictions in common between Rockhopper and APERO 

IGRs_summary_Rockhopper = IGRs_summary.dropna(subset=['Rockhopper_1'])

IGRs_summary_APERO = IGRs_summary.dropna(subset=['APERO_1'])


#%% save results




