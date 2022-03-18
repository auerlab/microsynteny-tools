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
    printf "Usage: $0 adjacent-genes max-nt gene-name [gene-name ...]\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# -lt 3 ]; then
    usage
fi
adjacent_genes=$1
max_nt=$2
shift
shift

make clean install > /dev/null
for species in Danio_rerio Oryzias_latipes Takifugu_rubripes; do
    for gene in $@; do
	gene_files="Regions/$species-$gene-*-$adjacent_genes-$max_nt.gff3"
	for file in $gene_files; do
	    if [ -e $file ]; then
		regen="$regen $file"
	    fi
	done
    done
done
for species in Mus_musculus Rattus_norvegicus Homo_sapiens; do
    for gene in $@; do
	gene_files="Regions/$species-$gene-*-$adjacent_genes-$max_nt.gff3"
	for file in $gene_files; do
	    if [ -e $file ]; then
		noregen="$noregen $file"
	    fi
	done
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
