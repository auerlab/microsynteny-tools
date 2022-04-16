#!/bin/sh -e

##########################################################################
#   Synopsis:
#       
#   Description:
#       Pipeline script for fish mammal comparison, running GFF
#       extraction, region intersects, and percent conserved analysis.
#
#   Arguments:
#       
#   Returns:
#
#   Examples:
#
#   Files:
#
#   Environment:
#
#   See also:
#       
#   History:
#   Date        Name        Modification
#   2022-04-04  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 adjacent-genes max-nt\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 2 ]; then
    usage
fi
adjacent_genes=$1
max_nt=$2

if [ ! -e Genes/de-tf-ortho.txt ] || [ ! -e Genes/non-de-tf-ortho.txt ]; then
    cd Genes/Orig
    printf "Extracting gene lists...\n"
    extract.sh
    cd ..
    printf "Adding biomart orthologs...\n"
    add-ortho.sh
fi

if [ -z $(ls Regions/*-$adjacent_genes-$max_nt.gff3) ]; then
    for list in de-tf-ortho.txt non-de-tf-ortho.txt; do
	gene-list-extract.sh --adjacent-genes $adjacent_genes Genes/$list
    done
    rm -f fish-mammal-percent.out # Must be regenerated if regions have changed
fi

outfile=fish-mammal-percent-$adjacent_genes.out
if [ ! -e $outfile ]; then
    fish-mammal-percent-changed.sh $adjacent_genes $max_nt \
	2>&1 | tee $outfile
else
    printf "fish-mammal-percent.out already exists.  Remove it to rerun the percent changed stage.\n"
fi
