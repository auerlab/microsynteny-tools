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
for species in Danio_rerio Oryzias_latipes Takifugu_rubripes; do
    file=Regions/$species-$gene-*.gff3
    regen="$regen $file"
done
for species in Mus_musculus Rattus_norvegicus Homo_sapiens; do
    file=Regions/$species-$gene-*.gff3
    noregen="$noregen $file"
done

if [ $(echo $regen | wc -w) -ge 2 ]; then
    printf "\n====================\n"
    printf "$gene\n"
    printf "====================\n\n"
    printf "Neighboring genes conserved among the original group:\n\n"
    ./ms-intersect $regen --diverged $noregen

    if [ $(echo $noregen | wc -w) -ge 2 ]; then
	printf "\nNeighboring genes conserved among the diverged group:\n\n"
	./ms-intersect $noregen
    fi
fi

