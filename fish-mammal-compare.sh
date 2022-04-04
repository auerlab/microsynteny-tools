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
    printf "Usage: $0\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 0 ]; then
    usage
fi
if [ ! -e Genes/de-tf-ortho.txt ] || [ ! -e Genes/non-de-tf-ortho.txt ]; then
    cd Genes/Orig
    printf "Extracting gene lists...\n"
    ./extract.sh
    cd ..
    printf "Adding biomart orthologs...\n"
    ./add-ortho.sh
fi

if [ ! -e Regions ]; then
    for list in de-tf-ortho.txt non-de-tf-ortho.txt; do
	./gene-list-extract.sh --adjacent-genes 2 Genes/$list
    done
else
    printf "Regions already exists.  Remove it to rerun the extract stage.\n"
fi

if [ ! -e Intersects ]; then
    for list in de-tf-ortho.txt non-de-tf-ortho.txt; do
	./gene-list-fish-mammal-intersect.sh 2 1000000 Genes/$list \
	    2>&1 | tee fish-mammal-intersect.out
    done
else
    printf "Intersects already exists.  Remove it to rerun the intersect stage.\n"
fi

if [ ! -e fish-mammal-percent.out ]; then
    for list in de-tf-ortho.txt non-de-tf-ortho.txt; do
	./fish-mammal-percent-changed.sh 2 1000000 \
	    2>&1 | tee fish-mammal-percent.out
    done
else
    printf "fish-mammal-percent.out already exists.  Remove it to rerun the percent changed stage.\n"
fi
