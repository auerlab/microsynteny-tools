#!/bin/sh -e

##########################################################################
#   Description:
#       Convenient wrapper around ms-stack.py, placing species in
#       evolutionary order for easier viewing.
#
#   Arguments:
#       gene-name (e.g. jun, gap43)
#       
#   History:
#   Date        Name        Modification
#   2022-02-12  Jason Bacon Begin
##########################################################################

usage()
{
    cat << EOM

Usage: $0 [--show-gene-lens] adjacent-genes max-nt \\
	gene-name [gene-name ...]

Example: $0 4 100000 jun

	 will display a stack of Regions/*-jun-*-4-1000000.gff3

EOM
    exit 1
}


##########################################################################
#   Main
##########################################################################

while [ $(echo $1 | cut -c 1,1) = '-' ]; do
    flags="$flags $1"
    shift
done

if [ $# -lt 3 ]; then
    usage
fi
adjacent_genes=$1
max_nt=$2
shift
shift
printf "\n$*\n\n"

for species in Danio_rerio Oryzias_latipes Takifugu_rubripes \
    Mus_musculus Rattus_norvegicus Homo_sapiens; do
    for gene in $@; do
	gene_files="Regions/$species-$gene-*-$adjacent_genes-$max_nt*.gff3"
	for file in $gene_files; do
	    if [ -e $file ]; then
		files="$files $file"
	    fi
	done
    done
done
./ms-stack.py $flags $files
