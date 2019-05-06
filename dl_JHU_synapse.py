# !/usr/bin/env python
### doanload xenograft fastqs from RNAseq folder syn13363874
### and also xenograft bams from syn11678418
### make sure you have a .synapseConfig file with login credentials

import synapseclient
import synapseutils
import sys
import re
import os
import glob
from synapseclient import Project, File, Folder, Activity
from synapseclient import Schema, Column, Table, Row, RowSet, as_table_columns

### login to synapse and get the flatstat file
syn = synapseclient.Synapse()
syn.login()
entity = syn.tableQuery("SELECT id, specimenID, name FROM syn13363852 WHERE assay = 'rnaSeq' AND fileFormat = 'fastq' AND transplantationType = 'xenograft'   ")
files = list(entity)

for i in files[0:len(files)]:
        print(i)
        synID = i[0]
        specimenID = i[1].split(" ")[0]
        filename = i[2]
        filename = filename[:20] +"~"+ filename[21:] 
        # print(filename)
        syn.get(synID, downloadLocation = "/home/xdoan/MPNST/01_dl_files/output")
        os.rename("/home/xdoan/MPNST/01_dl_files/output/"+filename, "/home/xdoan/MPNST/01_dl_files/output/"+specimenID+"_"+filename)
