#!/bin/sh -e

##########################################################################
#   Description:
#       Extract genes and their neighbors from a GFF, saving the
#       neighborhood to another GFF.
#
#   History:
#   Date        Name        Modification
#   2022-02-04  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 [file.gff3 ...]\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# -lt 1 ]; then
    files=$(ls GFF/*.gff3)
else
    files="$@"
fi

make clean all
mkdir -p Hoods
for gff in $files; do
    printf "\n========================================================\n"
    printf "$gff\n"
    printf "========================================================\n"
    time ./ms-extract --max-nt-distance 10000000 --output-dir Hoods $gff gap43
done