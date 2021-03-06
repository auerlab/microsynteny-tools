#!/usr/bin/env python3

##########################################################################
#   Synopsis:
#       feature-view file.gff3 [file.gff3 ...]
#
#   Description:
#       Generate plots for GFF3 files with a small number of features,
#       such as those generated by ms-shire.
#       
#   Arguments:
#       file.gff3   A GFF3 file containing a small number of features
#       
#   Returns:
#       0 on success, non-zero error codes otherwise
#
#   Examples:
#       feature-view danio_rerio-jun.gff3
#
#   See also:
#       ms-shire(1), align-species(1)
#       
#   History:
#   Date        Name        Modification
#   2022-02-03  Jason Bacon Begin
##########################################################################

import sys,os
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord

argc = len(sys.argv)
if ( argc < 2 ):
    print("Usage:", sys.argv[0], "file.gff3 [file.gff3 ...]\n")
    exit(1)

for i in range(1, argc):
    filename = sys.argv[i]
    features = []

    # Create list of features to plot
    x_min=9999999999999999
    x_max=0
    c = 0
    with open(filename) as infile:
        for line in infile:
            if line[0] != '#':
                cols = line.split("\t")
                start = int(cols[3])
                end = int(cols[4])
                if start < x_min: x_min = start
                if end > x_max: x_max = end
                strand = 1 if cols[6] == '+' else -1
                features.append(GraphicFeature(start=start, end=end,
                            strand=strand, color="#33bbbb", label=cols[1]))
    infile.close()
    print("Start coordinate =", x_min, "End coordinate = ", x_max)
    print(features)
    
    # Prepare for matplotlib
    record = GraphicRecord(sequence_length = x_max, features=features)
    ax, _ = record.plot(with_ruler = True, figure_width = 12,
                x_lim = [x_min - 1000, x_max + 1000],
                strand_in_label_threshold = 7)
    ax.set_title(filename)
    plt.tight_layout()  # Reduce margins
    plt.show()

    stem = os.path.splitext(filename)[0]
    png_filename = ".".join([stem,"png"])
    ax.figure.savefig(png_filename, bbox_inches='tight')
