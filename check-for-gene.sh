#!/bin/sh

##########################################################################
#   Description:
#       
#   History:
#   Date        Name        Modification
#   2022-03-01  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 gene-name GFF-file [GFF-file ...]\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# -lt 2 ]; then
    usage
fi
gene=$1
shift

for file in "$@"; do
    echo $file
    awk '$3 == "gene"' $file | grep -i "name=$gene"
done
