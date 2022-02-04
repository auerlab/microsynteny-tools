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

import sys,os

import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
from dna_features_viewer import BiopythonTranslator
from BCBio import GFF

argc = len(sys.argv)
if ( argc < 2 ):
    print("Usage:", sys.argv[0], "file.gff3 [file.gff3 ...]\n")
    exit(1)

# Can build a record manually or use translate_record() as below
#features=[
#    GraphicFeature(start=1000, end=1020, strand=+1, color="#ffd700",
#                   label="Small feature"),
#    GraphicFeature(start=1020, end=1500, strand=+1, color="#ffcccc",
#                   label="Gene 1 with a very long name"),
#    GraphicFeature(start=1400, end=1700, strand=-1, color="#cffccc",
#                   label="Gene 2"),
#    GraphicFeature(start=1600, end=1900, strand=+1, color="#ccccff",
#                   label="Gene 3")
#]
#record = GraphicRecord(sequence_length=2000, features=features)
#record.plot()
#plt.show()

for i in range(1, argc):
    filename = sys.argv[i]

    # Determine boundaries of the whole neighborhood
    x_min=9999999999999999
    x_max=0
    with open(filename) as infile:
        for line in infile:
            if line[0] != '#':
                cols = line.split("\t")
                start = int(cols[3])
                if start < x_min:
                    x_min = start
                end = int(cols[4])
                if end > x_max:
                    x_max = end
    infile.close()
    print("Start coordinate =", x_min, "End coordinate = ", x_max)
    
    # Convert GFF file to a matplotlib record
    record = BiopythonTranslator().translate_record(filename, filetype="gff")
    # 12 of what I have no idea, but it about spans the screen
    ax, _ = record.plot(figure_width=12, x_lim = [x_min, x_max], strand_in_label_threshold=7)
    ax.set_title(filename)
    plt.show()

    stem = os.path.splitext(filename)[0]
    png_filename = ".".join([stem,"png"])
    ax.figure.savefig(png_filename, bbox_inches='tight')
