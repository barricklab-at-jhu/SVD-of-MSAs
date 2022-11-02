# Singular Value Decomposition (SVD) for Protein Sequences & MSAs

The following repository is for the code associated with "Singular value decomposition of protein sequences as a method to visualize sequence and residue space" by A. Koenigs*, G. El Nesr*, and D. Barrick.

## Table of Contents

* [Overview](#overview)
* [Quickstart](#quickstart)
* [Requirements](#requirements)
* [Workflow](#workflow)
* [Generating Alignments](#alignments)
* [Singular Value Decomposition](#svd)
* [Webscraper](#webscraper)


<a name="overview"/>

## Overview

Singular value decomposition (SVD) of multiple sequence alignments (MSAs) is an important and rigorous way to extract information about the sequences and residues that make up those sequences, including evaluating consensus and covariance sequence features, identifying subgroups within a sequence family.  This information can be correlated to structure, function, and stability. The following repository contains scripts and notebooks that explain and perform SVD in a way that is intuitive and accessible to protein scientists. We then use k-means clustering to define possible clusters. These clusters are then defined by their gene ontology (GO) terms from UniProt to identify specific clusters.

<a name="quickstart"/>

## Quickstart

Clone the repo: ```git clone https://github.com/barricklab-at-jhu/SVD-of-MSAs.git```

<a name="requirements"/>

## Requirements

All of the code written runs on Python 3.6 or newer. Scripts were written using Python 3.6 on a MacOSX (Unix) system. 

To run the Jupyter notebooks, please find more information <a href='https://jupyter.org'>here</a>.


The following non-standard Python libraries are required to run the **remove_gaps code**:
* biopython

The following non-standard Python libraries are required for the **SVD code**:
* numpy
* matplotlib
* pandas
* seaborn
* scipy
* biopython
* sklearn
* opencv-python

The following non-standard Python libraries are required to run the **webscraping code**:
* pandas
* beautifulsoup4
* lxml

To install all of these packages, open a command-line interface and run 

```pip3 install <package-name>``` 

To learn more about pip, visit <a href='https://pip.pypa.io/en/stable/installation/'>here</a> to get started.

<a name = "workflow"/>

## Workflow

1. Generate an adequate Multiple Sequence Alignment (MSA)
2. Run the singular-value decomposition function. This can be accessed from the Jupyter Notebook.
3. Run the Web-scraping script for Uniprot.org.

You can find more in-depth information about each of the associated steps, their requirements and their workflow. 

<a name = "alignments"/>

## Generating Multiple Sequence Alignments (MSA)

Properly generating a multiple sequence alignment is critical for how SVD will function and perform. The following steps were taken to generate the multiple sequence alignments used in this repo.

1. **Download a sequence set from PFam.** These sequences contain UniProt identifiers which is used later in the webscraping portion of this project.
2. **Remove all gaps.** Some of these sequences contain gaps, so we remove them all. To run the script, open commandline, enter the directory of the python script, and run 
   ```python3 remove_gaps.py <sequences_file.txt>```
4. **Remove sequences with ≥90% identity.** To remove duplicate or extremely similar sequences, we use <a href="http://weizhong-lab.ucsd.edu/cdhit_suite/cgi-bin/index.cgi">cd-hit</a>. Note that the parameter for percent identity can be changed in the settings of cd-hit. 
5. **Align the sequence to generate a MSA file.** To align the sequences, we used <a href="https://mafft.cbrc.jp/alignment/software/">MAFFT</a>. However, there are various other tools out there that will also create a multiple sequence alignment such as <a href="http://www.drive5.com/muscle/">MUSCLE</a> and <a href="http://www.clustal.org/clustal2/">CLUSTAL</a>.

The resulting multiple sequence alignment should be in FASTA format. This sequence set is now ready for SVD.

<a name = "svd"/>

## Singular Value Decomposition 

The **Homeodomain** folder contains all of the associated files for the Homeodomain protein. To see where to start, check the README in that folder.

The **Ras** folder contains all of the associated files for the Ras protein. To see where to start, check the README in that folder.

<a name = "webscraper"/>

## Webscraper for Uniprot (webscraper.py)

The webscraper.py script will take in csv file with at least one column (with a column ID) and output an associated csv file. The file is expected to have columns that include the protein ID as found on <a href= "https://www.uniprot.org">Uniprot</a>. The script only handles one column at a time. 

### Steps to Running webscraper.py

1. Open a command line interface and enter the directory containing the webscrapper
2. Run the script by using the following command:

``` python3 webscraper.py <filepath_to_cluster_csv> <output_filename_csv> <column_ID>```

  For example, to run the webscraper on the "blue" column/cluster for Homeodomain, simply run 
  
``` python3 webscraper.py Homeodomain/HD_3_clusters.csv Homeodomain/HD_3_blue.csv blue```

### The Output

The output of the script will be a CSV file that contains all of the following information for each of the proteins, if available on the Uniprot website. Otherwise, there will be either a "nan" or blank cell in that spot. Removing "nan" is as simple as doing a Find & Replace for the term in Excel. For more information about each of these identifiers, check the Uniprot website.

* Protein ID
* Identified Length
* Gene ID
* Protein Name
* Organism Name
* Taxonomic ID
* Molecular Function - GO Annot
* Biological Processes - GO Annot
* Cellular Component - GO Annot
