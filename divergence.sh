#!/bin/sh -e

##########################################################################
#   Synopsis:
#       
#   Description:
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
#   2022-02-17  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 gene-name\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
gene=$1

make
./ms-divergence \
    Hoods/Danio_rerio-$gene.gff3 \
    Hoods/Oryzias_latipes-$gene.gff3 \
    Hoods/Takifugu_rubripes-$gene.gff3 \
    Hoods/Xenopus_tropicalis-$gene.gff3 \
    Hoods/Gallus_gallus-$gene.gff3 \
    Hoods/Mus_musculus-$gene.gff3 \
    Hoods/Rattus_norvegicus-$gene.gff3 \
    Hoods/Homo_sapiens-$gene.gff3
