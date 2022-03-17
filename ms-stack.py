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

import sys, os
from os import path
import matplotlib.pyplot as plt
import re

#############################################################################
#   Description:
#
#   Arguments:
#   
#   Returns:
#   
#   History: 
#   Date        Name        Modification
#   2022-03-15  Jason Bacon Begin

def usage():
    print("Usage: %s %s" % (sys.argv[0], "[--no-gene-lens] [--batch] file.gff3 [file.gff3 ...]"))
    sys.exit(1)

#############################################################################
#   Process command line args

show_gene_lens = batch_mode = True
if len(sys.argv) > 1:
    arg = 1
    while (arg < len(sys.argv)) and (sys.argv[arg][0] == '-'):
        if sys.argv[arg] == '--no-gene-lens':
            show_gene_lens = False
        elif sys.argv[arg] == '--batch':
            batch_mode = True
        arg += 1
else:
    usage()

if show_gene_lens:
    arrow_y_sep   = 12  # A little extra room between lines
else:
    arrow_y_sep = 10

# Plotting parameters for all files
arrow_len   = 100
arrow_x_sep = 20
text_height = 0.7

# Starting point for first file
files = 0
for filename in sys.argv[arg:]:
    if path.exists(filename):
        files += 1
if files < 1:
    usage()
arrow_y = (files - 1) * arrow_y_sep

# plt.show() makes the width too small, so genes overlap even though they
# are explicitly spaced out.  Can we fix just the width?
plt.rcParams["figure.figsize"] = (12, len(sys.argv) * 0.75)

print("%-18s %2s  %s\n" % ("Species", "Ch", "Genes"), end='')
last_goi = gois = ''
for filename in sys.argv[arg:]:
    basename = path.basename(filename)
    
    # This depends on filename format used by ms-extract.c
    c = re.split("[-.]", basename)
    species = c[0]
    goi = c[1]
    chrom = c[2]
    adjacent_genes = c[4]
    max_nt = c[5]
    
    if goi != last_goi:
        gois = gois + ' ' + goi
        last_goi = goi
    
    #############################################################################
    #   Parse file line by line

    if path.exists(filename):
        print("%-18s %2s " % (species, chrom), end='')
        with open(filename) as infile:
            # Determine orientation of GOI
            genes = 0
            for gff_line in infile:
                if gff_line[0] != '#':
                    genes += 1
                    cols = gff_line.split("\t")
                    gene = cols[1]
                    strand = cols[6]
                    if gene.lower() == goi.lower():
                        goi_strand = strand
            
            x_increment = arrow_len + arrow_x_sep
            if goi_strand == '+':
                arrow_left = 200
            else:
                arrow_left = 200 + (genes - 1) * x_increment
                x_increment = -x_increment
            arrow_right = arrow_left + arrow_len
                
            #   Plot genes
            infile.seek(0)
            genes = 0
            previous_end = 0
            for gff_line in infile:
                
                # pyplot doesn't respect text, so the species name will be
                # off the screen unless we plot an element that pyplot
                # cares about with low enough coordinates
                plt.plot([50, 50], [arrow_y, arrow_y])
                
                plt.text(0, arrow_y - text_height, species + goi_strand)
                # FIXME: Compute X position based on max number of neighbors
                plt.text(170, arrow_y - text_height, chrom)
                if gff_line[0] != '#':
                    cols = gff_line.split("\t")
                    gene = cols[1]
                    gene_start = int(cols[3])
                    gene_end = int(cols[4])
                    strand = cols[6]
                    trunc = gene[0:7:]
                    if len(gene) > len(trunc):
                        trunc = trunc + '*' # Indicate truncation
                    
                    if strand == goi_strand:
                        print(" %s+" % (gene), end='')
                        arrow_start = arrow_left
                        dx = arrow_len
                        text_offset = 6
                    else:
                        print(" -%s" % (gene), end='')
                        arrow_start = arrow_right
                        dx = -arrow_len
                        text_offset = 9

                    if goi_strand == '+':
                        ig_x = arrow_left - 20
                    else:
                        ig_x = arrow_left - 20 - x_increment
                    
                    if gene.lower() == goi.lower():
                        color='#DDDD11'
                    else:
                        color='#33bbbb'
                    
                    plt.arrow(arrow_start, arrow_y, dx, 0,
                              width=3, head_length=10,
                              length_includes_head=True,
                              head_width=4, facecolor=color, edgecolor='black');
                    plt.text(arrow_left + text_offset, arrow_y - text_height, trunc)
                    
                    # Gene length
                    if show_gene_lens:
                        gene_len = gene_end - gene_start
                        plt.text(arrow_left + text_offset, arrow_y + 3,
                                 str(int(gene_len / 100 + 5) / 10) + 'k')
                    
                    # Intergenic distance
                    if genes > 0:
                        gap = gene_start - previous_end
                        # Kb rounded to 1 decimal digit
                        plt.text(ig_x, arrow_y - 4.5,
                            str(int(gap / 100 + 5) / 10) + 'k')
                    
                    arrow_left += x_increment
                    arrow_right = arrow_left + arrow_len
                    previous_end = gene_end
                    genes += 1
        infile.close()
        print()

    arrow_y -= arrow_y_sep

plt.title(gois + ' neighborhoods and intergenic distances, genes = ' +
    adjacent_genes + ' max_nt = ' + max_nt)
plt.box(False)
plt.xticks([])
plt.yticks([])
os.makedirs("Stacks", exist_ok=True)
plt.savefig("Stacks/" + goi + '-' + adjacent_genes + '-' + max_nt + "-stack.png")
# Creates a new blank figure, so do after savefig()
if not batch_mode:
    plt.show()
