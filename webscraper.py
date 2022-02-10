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

#display list
def display_list(arr):
    return ";\n".join(arr)

#get data; handle error
def fetchdata(url):
    try:
        return urlopen(url)
    except:
        fetchdata(url)

#get UniProt Protein IDs
input_file = sys.argv[1]
output_file = sys.argv[2]
code = sys.argv[3]
df = pd.read_csv(input_file)
IDs = df[code]

#Get molecular functions
with open(output_file, 'w') as csvfile1:
    #Write header row for Protein_data file
    datawriter1 = csv.writer(csvfile1, delimiter=',')
    header1 =    [  "Protein ID", "Identified Length", "Gene ID", "Protein Name", "Organism Name",
                    "Taxonomic ID", "Molecular Function-GO Annot","Molecular Function-Keyword",
                    "Biological processes-GO Annot","Biological processes-Keywords",
                    "Cellular Component-Go Annot", "Cellular Component-Keywords",
                    "Disease-OMIM ID", "Disease-Keywords",
                    "Technical Terms-Keywords", "Polymorphism" ]
    datawriter1.writerow(header1)

    countdown = len(IDs)
    print('Starting scraping of ' + code + ' from ' + input_file)
    start_time = time.time()
    for i in IDs:
        countdown -= 1
        try:
            pid = i.split('/')[0]
            plen = int(i.split('/')[1].split('-')[1]) - int(i.split('/')[1].split('-')[0]) + 2
            #specify the url
            url = "https://www.uniprot.org/uniprot/" + str(pid)
            #Query the website
            page = fetchdata(url)
            #Parse the html, store it in Beautiful Soup format
            bsf = BeautifulSoup(page, "lxml")

            #get Gene ID
            gid = ""
            ext_data = bsf.find('table', class_='databaseTable GENOME')
            if(not (ext_data is None)):
                data = ext_data.findAll(text=True)
                if "GeneID" in data:
                    i = data.index("GeneID")
                    gid = data[i+2]

            #Molecular Function GO Annotation
            molecular_function_go = []
            ext_data = bsf.find('ul', class_='noNumbering molecular_function')
            if(not (ext_data is None)):
                for data in ext_data.findAll('li'):
                    cells = data.find(lambda tag: tag.name == "a" and (tag.has_attr("onclick")))
                    if(not(cells is None)):
                        cells = cells.find(text=True).strip()
                        molecular_function_go.append(cells)

            #Biological Processes GO Annotation
            biological_process_go = []
            ext_data = bsf.find('ul', class_='noNumbering biological_process')
            if(not (ext_data is None)):
                for data in ext_data.findAll('li'):
                    cells = data.find(lambda tag: tag.name == "a" and (tag.has_attr("onclick")))
                    if(not(cells is None)):
                        cells = cells.find(text=True).strip()
                        biological_process_go.append(cells)

            #Cellular Component GO Annotation
            cellular_component_go = []
            ext_data = bsf.find('div', id='table-go_annotation')
            if(not (ext_data is None)):
                lt = ext_data.find('ul', class_='noNumbering subcellLocations')
                for li in lt.findAll('li'):
                    cell = li.find('h6')
                    if(not(cell is None)):
                        cell_loc = cell.find(text=True)
                        cellular_component_go.append(cell_loc)

            #protein names and taxonomic identifiers
            organism_name = []
            tax_id = []
            org = 0
            ext_data = bsf.find('div', id="names_and_taxonomy")
            if(not (ext_data is None)):
                ext_data = ext_data.find('table')
                for row in ext_data:
                    data = row.findAll('td')
                    head = data[0].find(text=True)
                    vals = data[1].findAll('a')
                    val = [v.find(text=True) for v in vals]
                    v = list(filter(lambda x : x != ', ', vals))
                    val = []
                    tax = []
                    for kw in v:
                        kws = str(kw)
                        if kws[:18]=='<a href="/taxonomy' and org == 0:
                            val.append(re.sub('<[^<]+?>', '', kw.find(text=True)))
                            organism_name = val
                            org = org + 1
                        elif kws[:18]=='<a href="/taxonomy' and org == 1:
                            tax.append(re.sub('<[^<]+?>', '', kw.find(text=True)))
                            tax_id = tax
                            org = org + 1

            #protein names
            protein_name = 'Uncharacterized protein'
            ext_data = bsf.find('span', class_="recommended-name")
            if(not (ext_data is None)):
                for row in ext_data:
                    protein_name = row
                    break


            #cellular component keywords
            cellular_component_kw = []
            ext_data = bsf.find('div', class_='section ', id="subcellular_location")
            if(not (ext_data is None)):
                header = ext_data.find('table')
                if(not (header is None)):
                    head_data = header.findAll(text=True)
                    for h in head_data:
                        data = header.next_sibling
                        vals = data.findAll('a')
                        val = [v.find(text=True) for v in vals]
                        v = list(filter(lambda x : x != ', ', vals))
                        val = []
                        for kw in v:
                            kws = str(kw)
                            if kws[:18]=='<a href="/keywords':
                                val.append(kw.find(text=True))
                        cellular_component_kw = val

            #Keywords - Molecular Function and Biological processes
            molecular_function_kw = []
            biological_process_kw = []
            ext_data = bsf.find('table', class_='databaseTable')
            if(not (ext_data is None)):
                ext_data = ext_data.findAll('tr')
                for row in ext_data:
                    data = row.findAll('td')
                    head = data[0].find(text=True)
                    vals = data[1].findAll('a')
                    val = [v.find(text=True) for v in vals]
                    v = list(filter(lambda x : x != ', ', vals))
                    val = []
                    for kw in v:
                        kws = str(kw)
                        if kws[:18]=='<a href="/keywords':
                            val.append(kw.find(text=True))
                    if(head=="Molecular function"):
                        molecular_function_kw = val
                    if(head=="Biological process"):
                        biological_process_kw = val

            #Disease OMIM ID
            disease_omim_id = []
            ids = []
            ext_data = bsf.findAll('div', class_='diseaseAnnotation')
            if(not (ext_data is None)):
                for data in ext_data:
                    val = data.findAll('a')
                    for data1 in val:
                        data2 = data.findAll(text=True)
                        for j in data2:
                            if j[:8]=="See also":
                                ids.append(j[14:])
            ids = set(ids)
            disease_omim_id = list(ids)

            #Disease Keywords
            disease_kw = []
            ext_data = bsf.find('div', class_='section', id='pathology_and_biotech')
            if(not (ext_data is None)):
                heads = ext_data.findAll('h4')
                for head in heads:
                    data = head.findAll(text=True)
                    if 'Keywords - Disease' in data:
                        j = data.index('Keywords - Disease')
                        val = data[j].parent.parent
                        cells = val.next_sibling
                        vals = cells.findAll('a')
                        val = [v.find(text=True) for v in vals]
                        v = list(filter(lambda x : x != ', ', vals))
                        val = []
                        for kw in v:
                            kws = str(kw)
                            if kws[:18]=='<a href="/keywords':
                                val.append(kw.find(text=True))
                        disease_kw = val
                        break

            #Technical Terms - Keywords
            tech_term_kw = []
            ext_data = bsf.find('div', class_='section', id='miscellaneous')
            if(not (ext_data is None)):
                heads = ext_data.findAll('h4')
                for head in heads:
                    data = head.findAll(text=True)
                    if 'Keywords - Technical term' in data:
                        j = data.index('Keywords - Technical term')
                        val = data[j].parent.parent
                        cells = val.next_sibling
                        vals = cells.findAll('a')
                        val = [v.find(text=True) for v in vals]
                        v = list(filter(lambda x : x != ', ', vals))
                        val = []
                        for kw in v:
                            kws = str(kw)
                            if kws[:18]=='<a href="/keywords':
                                val.append(kw.find(text=True))
                        tech_term_kw = val
                        break

            #Polymorphism
            polymorphism = ""
            ext_data = bsf.find('div', class_='section', id='sequences')
            if(not (ext_data is None)):
                heads = ext_data.findAll('h4')
                for head in heads:
                    data = head.findAll(text=True)
                    if 'Polymorphism' in data:
                        j = data.index('Polymorphism')
                        val = data[j].parent.parent
                        polymorphism = val.next_sibling.find(text=True)
                        break
            #write data to Protein_data CSV file
            print(str(pid))
            datawriter1.writerow([pid, plen, gid, protein_name, display_list(organism_name), display_list(tax_id),
                                display_list(molecular_function_go), display_list(molecular_function_kw),
                                display_list(biological_process_go), display_list(biological_process_kw),
                                display_list(cellular_component_go), display_list(cellular_component_kw),
                                display_list(disease_omim_id), display_list(disease_kw),
                                display_list(tech_term_kw), polymorphism]) 
        except:
            pass

df1 = pd.read_csv(output_file)
df1['Protein Name'] = df1['Protein Name'].map(lambda x: x.replace('<strong>', '').replace('</strong>', ''))
for column in ['Molecular Function-GO Annot', 'Molecular Function-Keyword', 'Biological processes-GO Annot', 'Biological processes-Keywords', 'Cellular Component-Go Annot', 'Cellular Component-Keywords', 'Disease-OMIM ID', 'Disease-Keywords', 'Technical Terms-Keywords', 'Polymorphism']:
    df1[column] = df1[column].astype(str)
    df1[column] = df1[column].str.rsplit(';\n')
    df1[column] = df1[column].str.join(', ')
df1.to_csv(output_file, index=False)

print('SCRAPING OF ' + code.upper() + ' FROM ' + input_file.upper() + ' COMPLETED')
print("--- %s seconds ---" % (time.time() - start_time))
