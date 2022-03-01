#!/bin/sh -e

##########################################################################
#   Description:
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

make > /dev/null
for species in Danio_rerio Oryzias_latipes Takifugu_rubripes \
    Mus_musculus Rattus_norvegicus Homo_sapiens; do
    file=Hoods/$species-$gene.gff3
    list="$list $file"
done
if [ $(echo $list | wc -w) -ge 2 ]; then
    printf "\n$gene\n\n"
    ./ms-intersect $list
fi

printf '\n'
list=''
for species in Mus_musculus Rattus_norvegicus Homo_sapiens; do
    file=Hoods/$species-$gene.gff3
    list="$list $file"
done
if [ $(echo $list | wc -w) -ge 2 ]; then
    ./ms-intersect $list
fi
