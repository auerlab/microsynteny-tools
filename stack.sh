#!/bin/sh -e

##########################################################################
#   Description:
#       Convenient wrapper around ms-stack.py, placing species in
#       evolutionary order for easier viewing.
#
#   Arguments:
#       gene-name (e.g. jun, gap43)
#       
#   History:
#   Date        Name        Modification
#   2022-02-12  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 gene\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
gene=$1
printf "\n$gene\n\n"

./ms-stack.py \
    Regions/Danio_rerio-$gene-*.gff3 \
    Regions/Oryzias_latipes-$gene-*.gff3 \
    Regions/Takifugu_rubripes-$gene-*.gff3 \
    Regions/Xenopus_tropicalis-$gene-*.gff3 \
    Regions/Gallus_gallus-$gene-*.gff3 \
    Regions/Mus_musculus-$gene-*.gff3 \
    Regions/Rattus_norvegicus-$gene-*.gff3 \
    Regions/Homo_sapiens-$gene-*.gff3

