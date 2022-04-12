#!/bin/sh -e

##########################################################################
#   Description:
#       
#   History:
#   Date        Name        Modification
#   2022-04-12  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0\n"
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

if [ $# != 0 ]; then
    usage
fi

selection=''
while [ 0"$selection" != 0q ]; do
    clear
    cat << EOM

1.. Add latest Biomart orthologs (Genes/biomart-orthologs.tsv.xz) to lists
2.. Extract regions
3.. Visualize regions in scale
4.. View stacked regions across species (scale discarded)
5.. Find intersections between regions across fish and mammals
6.. Show percent changed stats for fish and mammals
    (Extract regions for de-tf-ortho.txt and non-de-tf-ortho.txt first)
Q.. Quit

EOM
    printf "Selection? "
    read selection
    case "$selection" in
    1)
	cd Genes
	./add-ortho.sh
	cd ..
	;;
    
    2)
	printf "Adjacent genes? [3] "
	read genes
	if [ 0"$genes" = 0 ]; then
	    genes=3
	fi
	printf "Max distance in NT? [1,000,000] "
	read max_nt
	if [ 0"$max_nt" = 0 ]; then
	    max_nt=1000000
	fi
	ls Genes/*.txt
	printf "Gene list? [Genes/de-tf-refined-ortho.txt] "
	read gene_list
	if [ 0"$gene_list" = 0 ]; then
	    gene_list="Genes/de-tf-refined-ortho.txt"
	fi
	./gene-list-extract.sh \
	    --adjacent-genes $genes --max-nt-distance $max_nt $gene_list
	;;
    
    3)
	printf "Gene name? "
	read gene
	ls Regions/*-$gene-*.gff3
	printf "Adjacent genes? "
	read adjacent
	printf "Max NT distance? "
	read max_nt
	./plot.sh $gene $adjacent $max_nt
	;;

    4)
	printf "Gene name? "
	read gene
	ls Regions/*-$gene-*.gff3
	printf "Adjacent genes? "
	read adjacent
	printf "Max NT distance? "
	read max_nt
	./stack.sh $adjacent $max_nt $gene
	;;
    
    5)
	printf "Adjacent genes? [3] "
	read genes
	if [ 0"$genes" = 0 ]; then
	    genes=3
	fi
	printf "Max distance in NT? [1,000,000] "
	read max_nt
	if [ 0"$max_nt" = 0 ]; then
	    max_nt=1000000
	fi
	ls Genes/*.txt
	printf "Gene list? [Genes/de-tf-refined-ortho.txt] "
	read gene_list
	if [ 0"$gene_list" = 0 ]; then
	    gene_list="Genes/de-tf-refined-ortho.txt"
	fi
	./gene-list-fish-mammal-intersect.sh $genes $max_nt $gene_list
	;;
    
    6)
	printf "Adjacent genes? [3] "
	read genes
	if [ 0"$genes" = 0 ]; then
	    genes=3
	fi
	printf "Max distance in NT? [1,000,000] "
	read max_nt
	if [ 0"$max_nt" = 0 ]; then
	    max_nt=1000000
	fi
	./fish-mammal-percent-changed.sh $genes $max_nt
	;;
    
    Q|q)
	exit 0
	;;
	
    *)
	printf "Invalid selection: $selection\n"
	;;
    esac
    pause
done
