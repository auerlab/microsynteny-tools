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
    printf "Usage: $0 gene-name [gene-name ...]\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# -lt 1 ]; then
    usage
fi

make > /dev/null
for species in Danio_rerio Oryzias_latipes Takifugu_rubripes; do
    for gene in $@; do
	file=Regions/$species-$gene-*.gff3
	if [ -e $file ]; then
	    regen="$regen $file"
	fi
    done
done
for species in Mus_musculus Rattus_norvegicus Homo_sapiens; do
    for gene in $@; do
	file=Regions/$species-$gene-*.gff3
	if [ -e $file ]; then
	    noregen="$noregen $file"
	fi
    done
done

if [ $(echo $regen | wc -w) -ge 2 ]; then
    printf "\n====================\n"
    printf "$*\n"
    printf "====================\n\n"
    printf "Neighboring genes conserved among the original group:\n\n"
    ./ms-intersect $regen --diverged $noregen

    if [ $(echo $noregen | wc -w) -ge 2 ]; then
	printf "\nNeighboring genes conserved among the diverged group:\n\n"
	./ms-intersect $noregen
    fi
fi
