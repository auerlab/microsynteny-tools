#!/bin/sh -e

##########################################################################
#   Description:
#       
#   History:
#   Date        Name        Modification
#   2022-03-07  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 gene-name adjacent-genes max-nt\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 3 ]; then
    usage
fi
gene=$1
adjacent=$2
max_nt=$3

for file in Regions/*-$gene-*-$adjacent-$max_nt.gff3; do
    ./ms-plot.py $file
done
