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

#############################################################################
#   Process command line args

if len(sys.argv) < 2:
    print("Usage: %s %s" % (sys.argv[0], "file.gff3 [file.gff3 ...]"))
    sys.exit(1)

unrep = []

for filename in sys.argv[1:]:
    basename = path.basename(filename)
    c = basename.split("-")
    
    #############################################################################
    #   Parse file line by line

    if path.exists(filename):
        print("%-20s" % c[0], end='')
        with open(filename) as infile:
            for line in infile:
                if line[0] != '#':
                    cols = line.split("\t")
                    gene = cols[1]
                    strand = cols[6]
                    if strand == '+':
                        print("%s+ " % gene, end='')
                    else:
                        print("-%s " % gene, end='')
        infile.close()
        print()
    else:
        unrep.append(filename)

if len(unrep) > 0:
    print("\nThe following files were not found:")
    for file in unrep:
        print(file)
