#!/bin/sh -e

##########################################################################
#   Description:
#       Show genes near ONE gene-of-interest (GOI) conserved for all fish,
#       conserved for all fish and each individual mammal, and conserved
#       for all mammals.
#
#       Regional GFF files must have been previously generated with
#       ms-extract [--adjecent-genes adjacent-genes] [--max-nt max-nt].
#
#       If gene names are not exactly the same between species, you can
#       specify orthologous GOIs as extra arguments after the primary GOI
#       in the 3rd argument.
#
#   Arguments:
#       adjacent-genes  Number of adjacent genes on each side of GOI
#       max-nt          Max distance in NT for adjacent genes
#       gene-name       GOI
#       alt-gene-name   Orthologous gene name
#       
#   History:
#   Date        Name        Modification
#   2022-02-17  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 adjacent-genes max-nt gene-name [alt-gene-name ...]\n"
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
    # $@ should contain alternate genes such as isl2b isl2
    for goi in $@; do
	# '-' is a separator, so don't allow it in gene names
	# ms-extract substitutes '_'
	goi=$(echo $goi | tr '-' '_')
	gene_files="Regions/$species-$goi-*-$adjacent_genes-$max_nt.gff3"
	for file in $gene_files; do
	    if [ -e $file ]; then
		regen="$regen $file"
	    fi
	done
    done
done

for species in Mus_musculus Rattus_norvegicus Homo_sapiens; do
    # $@ should contain alternate genes such as isl2b isl2
    for goi in $@; do
	# '-' is a separator, so don't allow it in gene names
	# ms-extract substitutes '_'
	goi=$(echo $goi | tr '-' '_')
	gene_files="Regions/$species-$goi-*-$adjacent_genes-$max_nt.gff3"
	for file in $gene_files; do
	    if [ -e $file ]; then
		noregen="$noregen $file"
	    fi
	done
    done
done

if [ $(echo $regen | wc -w) -ge 2 ]; then
    printf "\n========================================================\n"
    printf "$*\n"
    printf "========================================================\n\n"
    printf "Genes conserved among fish and mammals:\n\n"
    # echo $regen
    # echo $noregen
    ./ms-intersect Intersects-diverged $regen --diverged $noregen

    if [ $(echo $noregen | wc -w) -ge 2 ]; then
	printf "\nGenes conserved among mammals only:\n\n"
	# echo $noregen
	./ms-intersect Intersects-diverged $noregen
    fi
fi
