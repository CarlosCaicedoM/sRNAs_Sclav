#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: % Carlos Caicedo-Montoya
"""

from Bio import Entrez  
import pandas as pd  


def search_experiments_by_species(species):  
    Entrez.email = 'candres.caicedo@udea.edu.co'  # Make sure to put your email address here  
    handle = Entrez.esearch(db='sra', term=species, retmax=1000)  # retmax is the maximum number of results to retrieve  
    results = Entrez.read(handle)  
    handle.close()  
    return results['IdList']  

def get_experiment_information(experiment_id):  
    handle = Entrez.esummary(db='sra', id=experiment_id)  
    result = Entrez.read(handle)  
    handle.close()  
    return result  

def extract_data(xml_data, start_tag, end_tag): 
    start_pos = xml_data.find(start_tag) 
    if start_pos != -1: 
        # Mover la posición inicial al final de la etiqueta de apertura 
        start_pos += len(start_tag) 
        # Encontrar la posición de la etiqueta de cierre después de la etiqueta de apertura 
        end_pos = xml_data.find(end_tag, start_pos) 
        if end_pos != -1: 
            # Extraer el texto entre las etiquetas de inicio y final 
            data = xml_data[start_pos:end_pos] 
            return data 

    # Si no se encuentra alguna de las etiquetas, retornar None 
    return None 

def main():
    species = "Streptomyces clavuligerus" 
    results = search_experiments_by_species(species) 
    metada_species = [] 
    for id in results: 
        metadata = get_experiment_information(id) 
        xml_data = metadata[0]["ExpXml"] #Here it is most of the information 
        metada_species.append(xml_data) 
    
    data = { 
        'Title': [], 
        'Instrument': [], 
        'Total Spots': [], 
        'Total Bases': [], 
        'Total Size': [], 
        'Submitter Acc': [], 
        'Center Name': [], 
        'Experiment Acc': [], 
        'Taxid': [], 
        'Scientific Name': [], 
        'Library Strategy': [], 
        'Library Source': [], 
        'Library Selection': [], 
        'Library Layout': [], 
        'Bioproject': [], 
        'Biosample': [] } 
    
    for xml_data in metada_species: 
        # Extract el título 
        title = extract_data(xml_data, '<Title>', '</Title>') 
        #Extract Instrument 
        instrument = extract_data(xml_data, '/><Instrument', '/><L') 
        #Extract total_spots 
        total_spots = extract_data(xml_data, 'total_spots="', '" total_bases') 
        #Extract total_bases 
        total_bases = extract_data(xml_data, 'total_bases="', '" total_size') 
        #Extract total_size 
        total_size = extract_data(xml_data, 'total_size="', '" load_done') 
        #Extract Submitter acc 
        submitter_acc = extract_data(xml_data, 'Submitter acc="', '" center_name=') 
        #Extract center_name 
        center_name = extract_data(xml_data, 'center_name="', '" contact_name') 
        #Extract Experiment acc 
        experiment_acc = extract_data(xml_data, 'Experiment acc="', '" ver=') 
        #Extract taxid 
        taxid = extract_data(xml_data, 'taxid="', '" ScientificName=') 
        #Extract ScientificName 
        ScientificName = extract_data(xml_data, 'ScientificName="', '"/><Sample') 
        #Extract LIBRARY_STRATEGY 
        LIBRARY_STRATEGY = extract_data(xml_data, 'LIBRARY_STRATEGY>', '</LIBRARY_STRATEGY>') 
        #Extract LIBRARY_SOURCE 
        LIBRARY_SOURCE = extract_data(xml_data, 'LIBRARY_SOURCE>', '</LIBRARY_SOURCE>') 
        #Extract LIBRARY_SELECTION 
        LIBRARY_SELECTION = extract_data(xml_data, 'LIBRARY_SELECTION>', '</LIBRARY_SELECTION') 
        #Extract LIBRARY_LAYOUT 
        LIBRARY_LAYOUT = extract_data(xml_data, 'LIBRARY_LAYOUT> <', '/> </LIBRARY_LAYOUT') 
        #Extract Bioproject 
        bioproject = extract_data(xml_data, 'Bioproject>', '</Bioproject') 
        #Extract Biosample 
        biosample = extract_data(xml_data, 'Biosample>', '</Biosample') 
      
        # Agregar los datos al diccionario 
        data['Title'].append(title) 
        data['Instrument'].append(instrument) 
        data['Total Spots'].append(total_spots) 
        data['Total Bases'].append(total_bases) 
        data['Total Size'].append(total_size) 
        data['Submitter Acc'].append(submitter_acc) 
        data['Center Name'].append(center_name) 
        data['Experiment Acc'].append(experiment_acc) 
        data['Taxid'].append(taxid) 
        data['Scientific Name'].append(ScientificName) 
        data['Library Strategy'].append(LIBRARY_STRATEGY) 
        data['Library Source'].append(LIBRARY_SOURCE) 
        data['Library Selection'].append(LIBRARY_SELECTION) 
        data['Library Layout'].append(LIBRARY_LAYOUT) 
        data['Bioproject'].append(bioproject) 
        data['Biosample'].append(biosample) 
    
    
    run_acc = [] 
    for id in results: 
        metadata = get_experiment_information(id) 
        run_data = metadata[0]["Runs"] #Here it is most of the information 
        run_acc.append(run_data) 
    
    run_acc_list = []
    for i in run_acc: 
        # Extract run acc 
        run = extract_data(i, 'Run acc="', '" total_spots') 
        run_acc_list.append(run)
    
    data['run_acc'] = run_acc_list
    
    df = pd.DataFrame(data)
    transcriptomics = df.loc[df['Library Source'] == 'TRANSCRIPTOMIC']
    rna_seq = transcriptomics.loc[transcriptomics['Library Strategy'] == 'RNA-Seq']
    rna_seq.to_csv("experiment_data.csv", index=False)

if __name__ == "__main__": 
    main() 