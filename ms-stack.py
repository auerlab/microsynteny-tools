#!/usr/bin/env python3

##########################################################################
#   Synopsis:
#       align-species file.gff3 [file.gff3 ...]
#
#   Description:
#       Display the gene neighborhoods of multiple files on top of each
#       other.  This is most useful for comparing the gene neighbohood
#       of a particular gene across two or more species.
#
#       The GFF3 files are generally output from ms-shire, containing
#       a small set of consecutive gene features surrounding a gene of
#       interest.
#
#   Arguments:
#       file.gff3   A GFF3 file containing a gene neighborhood
#       
#   Returns:
#       0 on success, non-zero error codes if a failure occurs
#
#   Examples:
#       align-species Regions/Danio_rerio-jun.gff3 \\\\
#           Regions/Takifugu_rubripes-jun.gff3 Regions/Xenopus_tropicalis-jun.gff3 \\\\
#           Regions/Mus_musculus-jun.gff3
#       Danio   si:dkey-239i20.4 plpp6 prdx6 jun caiap si:dkey-86e18.1 faslg
#       Fugu    unnamed plpp6 prdx6 jun si:dkey-86e18.1 faslg fam20b
#       Xenopus mysm1 unnamed pan3 jun fggy XB5864909 [provisional:rnf170] hook1
#       Mouse   Tek Eqtn Mysm1 Jun Fggy Hook1 Cyp2j13
#
#   See also:
#       ms-shire(1), feature-view(1)
#       
#   History:
#   Date        Name        Modification
#   2022-02-06  Jason Bacon Begin
##########################################################################

import sys
from os import path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

#############################################################################
#   Process command line args

if len(sys.argv) < 2:
    print("Usage: %s %s" % (sys.argv[0], "file.gff3 [file.gff3 ...]"))
    sys.exit(1)

unrep = []

bar_len = 60
bar_sep = 20
bar_y   = 0

# plt.show() makes the width too small, so genes overlap even though they
# are explicitly spaced out.  Can we fix just the width?
plt.rcParams["figure.figsize"] = (14, 5)

print("%-18s %2s  %s\n" % ("Species", "Ch", "Genes"), end='')
for filename in sys.argv[1:]:
    basename = path.basename(filename)
    
    # This depends on filename format used by ms-extract.c
    c = basename.split("-")
    species = c[0]
    goi = c[1]
    chrom = c[2]
    
    #############################################################################
    #   Parse file line by line

    if path.exists(filename):
        bar_left = 120
        bar_right = bar_left + bar_len
        print("%-18s %2s " % (species, chrom), end='')
        with open(filename) as infile:
            for gff_line in infile:
                plt.text(0, bar_y - 1, species)
                plt.text(96, bar_y - 1, chrom)
                if gff_line[0] != '#':
                    cols = gff_line.split("\t")
                    gene = cols[1]
                    start = cols[4]
                    strand = cols[6]
                    if strand == '+':
                        print(" %s+" % (gene), end='')
                        trunc = gene[0:8:] + ' +'
                    else:
                        print(" -%s" % (gene), end='')
                        trunc = '- ' + gene[0:8:]
                    plt.plot([bar_left, bar_right], [bar_y, bar_y],
                             linestyle='-', linewidth=16, color='#66CCCC')
                    plt.text(bar_left + 4, bar_y - 1, trunc, color='black')
                    plt.text(bar_left + 2, bar_y - 6,
                             str(int(int(start) / 1000)) + 'k')
                    bar_left = bar_right + bar_sep
                    bar_right = bar_left + bar_len
        infile.close()
        print()
    else:
        unrep.append(filename)

    bar_y += 10

plt.title(goi + ' neighborhoods')
plt.box(False)
plt.xticks([])
plt.yticks([])
plt.show()

#if len(unrep) > 0:
#    print("\nThe following files were not found:\n")
#    for file in unrep:
#        print(file)

