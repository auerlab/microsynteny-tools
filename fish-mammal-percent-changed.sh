#!/bin/sh -e

##########################################################################
#   Description:
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
    printf "Usage: $0 adjacent-genes max-nt\n"
    exit 1
}


##########################################################################
#   Function description:
#       Pause until user presses return
##########################################################################

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}

##########################################################################
#   Main
##########################################################################

if [ $# -lt 2 ]; then
    usage
fi
adjacent_genes=$1
max_nt=$2

make clean install > /dev/null

regen_species="Danio_rerio Oryzias_latipes Takifugu_rubripes"
noregen_species="Mus_musculus Rattus_norvegicus Homo_sapiens"

for gene_list in Genes/de-tf.txt Genes/non-de-tf.txt; do
    rm -f counts.txt
    for goi in $(cat $gene_list); do
	orthologs=$(echo $goi | tr '|' ' ')
	regen=''
	for species in $regen_species; do
	    for gene in $orthologs; do
		gene_files="Regions/$species-$gene-*-$adjacent_genes-$max_nt.gff3"
		for file in $gene_files; do
		    if [ -e $file ]; then
			regen="$regen $file"
		    fi
		done
	    done
	done
	
	if [ $(echo $regen | wc -w) -ge 2 ]; then
	    ./ms-intersect $regen > temp.out
	    awk '$1 == "Genes:" { print $2, $4, $6 }' temp.out >> counts.txt
	fi
    done

    printf "\nOverall conservation rate for $gene_list, $adjacent_genes, $max_nt\n"
    printf "$regen_species:\n"
    awk -f conserved-stats.awk counts.txt
    
    rm -f counts.txt
    for goi in $(cat $gene_list); do
	orthologs=$(echo $goi | tr '|' ' ')
	noregen=''
	for species in $noregen_species; do
	    for gene in $orthologs; do
		gene_files="Regions/$species-$gene-*-$adjacent_genes-$max_nt.gff3"
		for file in $gene_files; do
		    if [ -e $file ]; then
			noregen="$noregen $file"
		    fi
		done
	    done
	done
    
	if [ $(echo $noregen | wc -w) -ge 2 ]; then
	    ./ms-intersect $noregen > temp.out
	    awk '$1 == "Genes:" { print $2, $4, $6 }' temp.out >> counts.txt
	fi
    done

    printf "\nOverall conservation rate for $gene_list, $adjacent_genes, $max_nt\n"
    printf "$noregen_species:\n"
    awk -f conserved-stats.awk counts.txt
done
rm -f temp.out counts.txt
