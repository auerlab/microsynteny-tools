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

case $# in
1)
    gene_list=$1
    files="$(ls GFF/*.gff3)"
    ;;
2)
    gene_list=$1
    shift
    files="$@"
    ;;
*)
    usage
    ;;
esac

make clean all
awk '{ print $2 }' DETF_refined.tsv > DETF-refined.txt
mkdir -p Hoods
for gff in $files; do
    printf "\n========================================================\n"
    printf "$gff\n"
    printf "========================================================\n"
    time ./ms-extract --max-nt-distance 10000000 --output-dir Hoods \
	$gff $gene_list
done
