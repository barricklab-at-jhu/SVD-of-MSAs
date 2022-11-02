# Modified source code from https://github.com/Aksh77/Bio-Scraper/blob/master/UniProt-Scraper/scraper.py
import os
import re
import sys
import csv
import pandas as pd
import pprint
import time
import urllib
from urllib.request import urlopen
from bs4 import BeautifulSoup
import requests

# Display list
def display_list(arr):
    return ";\n".join(arr)

# Handle
def fetchdata(url):
    try:
        return urlopen(url)
    except:
        fetchdata(url)

# Get UniProt Protein IDs
input_file = sys.argv[1]
output_file = sys.argv[2]
code = sys.argv[3]
df = pd.read_csv(input_file)
IDs = df[code]

with open(output_file, 'w') as csvfile1:
    # Write header row for Protein_data file
    datawriter1 = csv.writer(csvfile1, delimiter=',')
    header1 =    [  "Protein ID", "Identified Length", "Gene ID", "Protein Name", "Organism Name",
                    "Taxonomic ID", "Molecular Function-GO Annot",
                    "Biological processes-GO Annot",
                    "Cellular Component-Go Annot"]
    datawriter1.writerow(header1)

    countdown = len(IDs)
    print('Starting scraping of ' + code + ' from ' + input_file)
    start_time = time.time()
    
    for i in IDs[:10]:
        countdown -= 1
        try:
            pid = i.split('/')[0]
            plen = int(i.split('/')[1].split('-')[1]) - int(i.split('/')[1].split('-')[0]) + 2
            
            # Specify URL
            url = "https://www.uniprot.org/uniprotkb/" + str(pid) + ".xml"
            req = requests.get(url)

            # Query Website
            bsf = BeautifulSoup(req.content, 'xml')
            
            # Gene ID
            gid = ""
            ext_data = bsf.find("name", type="ORF")
            if ext_data is not None:
                gid = ext_data.string

            # Molecular Function GO Annotation
            # Cellular Component GO Annotation
            # Biological Processes GO Annotation
            molecular_function_go = []
            cellular_component_go = []
            biological_process_go = []
            
            ext_data = bsf.find_all('dbReference')
            if ext_data is not None:
                for ed in ext_data:
                    if ed.find('property', value="ECO:0007669") is not None:
                        term = ed.find('property', type="term")["value"]
                        if term[0] == 'F':
                            molecular_function_go.append(term[2:])
                        elif term[0] == 'C':
                            cellular_component_go.append(term[2:])
                        elif term[0] == 'P':
                            biological_process_go.append(term[2:])

            # Taxonomic Identifier
            ext_data = bsf.find_all('dbReference', type="NCBI Taxonomy")
            tax_id = ext_data[0]['id']

            # Organism
            organism_name = ''
            ext_data = bsf.find('name', type="scientific")
            if ext_data is not None:
                organism_name = ext_data.string
            
            ext_data = bsf.find_all('name', type="common")
            if ext_data is not None:
                for ed in ext_data:
                    organism_name += ' ('+ ed.string +')'
            
            ext_data = bsf.find_all('name', type="synonym")
            if ext_data is not None:
                for ed in ext_data:
                    organism_name += ' ('+ ed.string +')'

            # Protein Name
            protein_name = 'Uncharacterized protein'
            ext_data = bsf.find('recommendedName')
            if ext_data is not None:
                ed = ext_data.find('fullName')
                if ed is not None:
                    protein_name = ed.string

            #write data to Protein_data CSV file
            datawriter1.writerow([
                pid, 
                plen, 
                gid, 
                protein_name, 
                organism_name, 
                tax_id,
                display_list(molecular_function_go),
                display_list(biological_process_go),
                display_list(cellular_component_go)]) 

        except:
            pass

df1 = pd.read_csv(output_file)
df1['Protein Name'] = df1['Protein Name'].map(lambda x: x.replace('<strong>', '').replace('</strong>', ''))
for column in ['Molecular Function-GO Annot', 'Biological processes-GO Annot', 'Cellular Component-Go Annot']:
    df1[column] = df1[column].astype(str)
    df1[column] = df1[column].str.rsplit(';\n')
    df1[column] = df1[column].str.join(', ')
df1.to_csv(output_file, index=False)

print('SCRAPING OF ' + code.upper() + ' FROM ' + input_file.upper() + ' COMPLETED')
print("--- %s seconds ---" % (time.time() - start_time))
