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
    printf "Usage: $0 gene-list.txt\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
gene_file="$1"

for gene in $(cat $gene_file); do
    ./stack.sh $(echo $gene | tr '|' ' ')
done
