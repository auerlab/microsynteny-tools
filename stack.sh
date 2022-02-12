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

./ms-stack.py \
    Hoods/Danio_rerio-$gene.gff3 \
    Hoods/Takifugu_rubripes-$gene.gff3 \
    Hoods/Xenopus_tropicalis-$gene.gff3 \
    Hoods/Mus_musculus-$gene.gff3

