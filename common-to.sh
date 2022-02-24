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
    Takifugu_rubripes \
    Oryzias_latipes; do
    file=Hoods/$species-$gene.gff3
    if [ -e $file ]; then
	list="$list $file"
    fi
done
make
if [ $(echo $list | wc -w) -ge 2 ]; then
    printf "\n$gene\n\n"
    ./ms-common-to $list
fi
