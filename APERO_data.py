#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: % Carlos Caicedo-Montoya
"""

import os  
import pandas as pd  

 
 
def process_files(file1, file2): 
# Leer el primer archivo CSV y filtrar según los criterios dados 
    df1 = pd.read_csv(os.path.join(apero_folder, file1)) 
    df1_coverage = df1[df1['freq'] >= 10] 
    df1_internal = df1_coverage[df1_coverage['Class'] != 'Ai'] 
        
    # Leer el segundo archivo CSV y filtrar según los criterios dados 
    df2 = pd.read_csv(os.path.join(apero_folder, file2)) 
    df2_coverage = df2[df2['freq'] >= 10] 
    df2_internal = df2_coverage[df2_coverage['Class'] != 'Ai']
    
    #Los small RNAs son menores que 500 y además en este caso mi pipeline  
    #para detectar sRNAs solo considera aquellos menores que esta distancia 
    df1_lg = df1_internal[df1_internal['lg'] <= 500]  
    df2_lg = df2_internal[df2_internal['lg'] <= 500]  
        
    return df1_lg, df2_lg 

def merge_rows(df_internal):
    # Sort the DataFrame by the 'Position' column
    df_sorted = df_internal.sort_values(by='Position').reset_index(drop=True)

    # Initialize a list to store merged rows
    merged_rows = []

    # Initialize variables to keep track of the merging process
    start_index = 0
    current_index = 1
    max_lg_row = df_sorted.iloc[start_index]

    # Iterate over all rows in the DataFrame
    while current_index < len(df_sorted):
        # Calculate the difference in the 'Position' column between the current and next element
        position_diff = df_sorted.loc[current_index, 'Position'] - df_sorted.loc[current_index - 1, 'Position']

        # Check if the difference is less than 10
        if position_diff < 10:
            # Update the row with the maximum 'lg' if necessary
            if df_sorted.loc[current_index, 'lg'] > max_lg_row['lg']:
                max_lg_row = df_sorted.iloc[current_index]
        else:
            # Append the row with the maximum 'lg' to the merged list
            merged_rows.append(max_lg_row)

            # Update the start index and reset the row with the maximum 'lg'
            start_index = current_index
            max_lg_row = df_sorted.iloc[start_index]

        # Increment the current index
        current_index += 1

    # Add the last row with the maximum 'lg' to the merged list
    merged_rows.append(max_lg_row)

    # Create a new DataFrame with the merged rows
    result = pd.DataFrame(merged_rows)

    return result

def merge_results(results1, results2):
    # Initialize the DataFrame to store the merged results
    merged_results = pd.DataFrame()

    # Iterate over each row of the first DataFrame (results1)
    for index, row1 in results1.iterrows():
        position1 = row1['Position']
        lg1 = row1['lg']

        # Filter corresponding rows in the second DataFrame (results2)
        matching_rows = results2[
            (results2['Position'] >= (position1 - 10)) &
            (results2['Position'] <= (position1 + 10))
        ]

        # Find the row in matching_rows with the maximum "lg" value
        if not matching_rows.empty:
            max_lg_row = matching_rows.loc[matching_rows['lg'].idxmax()]

            # Compare the "lg" values of row1 and max_lg_row
            if lg1 > max_lg_row['lg']:
                merged_results = merged_results.append(row1)
            else:
                merged_results = merged_results.append(max_lg_row)

    # Reset the indices of the merged DataFrame
    merged_results.reset_index(drop=True, inplace=True)

    return merged_results


my_experiment = "S-clav_RNA-seq_sRNas_OMEGA" 
cur_dir = os.getcwd() 
apero_folder= os.path.join(cur_dir, "results_rna_seq", my_experiment, "APERO")  

#Open results from replicate 1  
replicates1 = [file for file in os.listdir(apero_folder) if (file.endswith("-1.csv")) ] 
replicates1 = set(replicates1) 
sorted_replicates1 = sorted(replicates1) 
 #Open results from replicate 2 
replicates2 = [file for file in os.listdir(apero_folder) if (file.endswith("-2.csv")) ] 
replicates2 = set(replicates2) 
sorted_replicates2 = sorted(replicates2) 

# Inicializar una lista para almacenar los DataFrames procesados 
processed_dataframes = [] 
# Iterar sobre cada par de nombres de archivo y procesarlos 
for file1, file2 in zip(sorted_replicates1, sorted_replicates2): 
    df1_internal, df2_internal = process_files(file1, file2) 
    processed_dataframes.append((df1_internal, df2_internal)) 
 
dataframes_filtered_size = [element for element in processed_dataframes if all(len(df) > 0 for df in element)]

#Delete consecutive rows with the sama information
merged_rows = []
for i in dataframes_filtered_size:
    # Aplicar la función a cada DataFrame en el elemento y guardar el resultado en una nueva lista
    result_i = [merge_rows(df) for df in i]
    # Agregar el resultado a la nueva lista
    merged_rows.append(result_i)

# merge the rows of the replicates that share the same annotations
merged_replicates = []
for i in merged_rows:
    merged_replicates.append(merge_results(i[0], i[1]))


results_concatenated = pd.concat(merged_replicates, ignore_index=True)
final_result=merge_rows(results_concatenated)
