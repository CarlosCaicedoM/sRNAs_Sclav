#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 11:41:56 2024

@author: usuario
"""
   

    
#%% Re-do the analysis considering the genes flanking the IGR
import pandas as pd
import os

bioinformatics_predictions = pd.read_table("results_sRNA_clade2_pangenome_CDS/reference_results/putative_sRNAs_rfam_blast_infernal_CP.csv", sep=",")
bioinformatics_predictions = bioinformatics_predictions.rename(columns={bioinformatics_predictions.columns[0]: 'header'})


rockhopper_results_CHR = pd.read_table("results_rna_seq/S-clavuligerus_RNA-Seq-data_SRA/Rockhopper3/NZ_CP027858_transcripts.txt", sep="\t") 

desired_columns = ['Transcription Start', 'Translation Start', 'Translation Stop',
                   'Transcription Stop', 'Strand', 'Name', 'Synonym', 'Product']

# Filter columns that start with 'Raw'
raw_columns = [col for col in rockhopper_results_CHR.columns if col.startswith('Raw')]

# Combine the lists of desired columns and columns starting with 'Raw'
columns_to_include = desired_columns + raw_columns

# Create a new DataFrame with the specified columns
rockhopper_results_CHR_filtered = rockhopper_results_CHR[columns_to_include]


rockhopper_results_CHR_filtered_nonAS = rockhopper_results_CHR_filtered[~rockhopper_results_CHR_filtered['Product'].str.startswith('antisense:')]
df = rockhopper_results_CHR_filtered_nonAS
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
filtered_df = IGR_sRNAS_rockhopper[IGR_sRNAS_rockhopper[raw_columns2].gt(0).sum(axis=1) > 10]
final_df = filtered_df[(filtered_df[raw_columns2] > 20).sum(axis=1) >= 20]
# Counter for numbering "Predicted RNA" entries
counter = 1

# Iterate over the rows of the DataFrame
for index, row in final_df.iterrows():
    # Check if the entry in the "Synonym" column is "Predicted RNA"
    if row['Synonym'] == 'predicted RNA':
        # Generate the new value for the "Synonym" entry
        new_value = f'predicted_RNA_{counter}'
        # Update the entry in the DataFrame
        final_df.at[index, 'Synonym'] = new_value
        # Increment the counter for the next entry
        counter += 1    


#%% Compare bioinformatics predictions and Rockhopper

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
    for j in final_df.itertuples():
        if (flanking_genes == j.Flanking):
            conserved_region = (j.Transcription_Start, j.Transcription_Stop, j.Strand, j.Synonym)
            predictions_in_Rockhopper[key].append(conserved_region)

 
df2_predictions_in_Rockhopper= pd.DataFrame.from_dict(predictions_in_Rockhopper, orient="index")

counter = 0

# Iterate over the dictionary and count elements that are not empty lists
for value_list in predictions_in_Rockhopper.values():
    if value_list:  # Check if the list has elements (i.e., not empty)
        counter += 1

verified_sRNAs  =  counter

print(verified_sRNAs)

#Concatenate bioinformatics_predictions with Rockhopper Predictions
bioinformatics_predictions.set_index('header', inplace=True)
putative_sRNAs_rfam_blast_infernal_CP_rockhopper = pd.concat([bioinformatics_predictions, 
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


#Concatenate bioinformatics_predictions with Rockhopper Predictions
putative_sRNAs_rfam_blast_infernal_CP_rockhopper_apero = pd.concat([putative_sRNAs_rfam_blast_infernal_CP_rockhopper, 
                                                df2_predictions_in_APERO], axis = 1)






#comparar apero vs rockhopper




# crear una sola tabla donde se vean si la region fue predicha por
# apero, rockhopper o mis predicciones basado en las pedicciones
#como por ejemplo decir, validado por...


#Restringir predicciones de rocjhopper para tener solo aquellas que aparezcan en 
#al menos 10 predicciones experimentales




