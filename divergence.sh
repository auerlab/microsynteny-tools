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

for species in \
    Danio_rerio \
    Oryzias_latipes \
    Takifugu_rubripes \
    Xenopus_tropicalis \
    Gallus_gallus \
    Mus_musculus \
    Rattus_norvegicus \
    Homo_sapiens; do
    file=Regions/$species-$gene.gff3
    if [ -e $file ]; then
	list="$list $file"
    fi
done
make
if [ $(echo $list | wc -w) -ge 2 ]; then
    printf "\n$gene\n\n"
    ./ms-divergence $list
fi
