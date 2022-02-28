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
    printf "Usage: $0 gene-list.txt [file.gff3 ...]\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# = 0 ]; then
    usage
fi

gene_list=$1
shift

if [ $# != 0 ]; then
    while [ 0$(echo $1 | cut -c 1) = 0'-' ]; do
	# All flags have one argument (for now)
	flags="$flags $1 $2"
	shift
	shift
    done
fi

if [ $# = 0 ]; then
    files=$(ls GFF/*.gff3)
else
    files="$@"
fi
printf "%s\n" "$flags"
printf "$files\n"

make clean all
mkdir -p Hoods
for gff in $files; do
    printf "\n========================================================\n"
    printf "$gff\n"
    printf "========================================================\n"
    time ./ms-extract --output-dir Hoods $flags $gff $gene_list
done
