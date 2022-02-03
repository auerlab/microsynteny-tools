#!/usr/bin/env python3

#############################################################################
#   Main Program
#
#   Description:
#
#   Arguments:
#   
#   Returns:
#   
#   History: 
#   Date        Name        Modification
#   2022-02-03  Jason Bacon Begin

# FIXME: biopython seems to require py-sqlite3
# FIXME: Docs don't show explicit plt.show()
# FIXME: ImportError("Please install the bcbio-gff library to parse GFF data")

# Quick-start:
# https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/

# import sys,os

import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
from dna_features_viewer import BiopythonTranslator
from BCBio import GFF

# Can build a record manually or use translate_record() as below
features=[
    GraphicFeature(start=1000, end=1020, strand=+1, color="#ffd700",
                   label="Small feature"),
    GraphicFeature(start=1020, end=1500, strand=+1, color="#ffcccc",
                   label="Gene 1 with a very long name"),
    GraphicFeature(start=1400, end=1700, strand=-1, color="#cffccc",
                   label="Gene 2"),
    GraphicFeature(start=1600, end=1900, strand=+1, color="#ccccff",
                   label="Gene 3")
]
#record = GraphicRecord(sequence_length=2000, features=features)
#record.plot()
#plt.show()

file="jun-hood.gff3"
handle = open(file)
for rec in GFF.parse(handle):
    print(rec)
handle.close()

record = BiopythonTranslator().translate_record("jun-hood.gff3", filetype="gff")
print(record)
ax, _ = record.plot(x_lim = [15550000, 16000000], strand_in_label_threshold=7)
plt.show()
