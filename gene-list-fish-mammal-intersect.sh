#!/bin/sh -e

##########################################################################
#   Description:
#       
#   History:
#   Date        Name        Modification
#   2022-03     Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 adjacent-genes max-nt gene-list.txt\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 3 ]; then
    usage
fi
adjacent_genes=$1
max_nt=$2
gene_file="$3"

for gene in $(cat $gene_file); do
    # Each line in gene_file may contain alternates such as isl2b|isl2
    ./fish-mammal-intersect.sh $adjacent_genes $max_nt $(echo $gene | tr '|' ' ')
done
