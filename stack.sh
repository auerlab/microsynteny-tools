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
    printf "Usage: $0 gene-name [gene-name ...]\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# -lt 1 ]; then
    usage
fi
printf "\n$*\n\n"

for species in Danio_rerio Orizias_latipes Takifugu_rubripes \
    Mus_musculus Rattus_norvegicus Homo_sapiens; do
    for gene in $@; do
	file=Regions/$species-$gene-*.gff3
	if [ -e $file ]; then
	    files="$files $file"
	fi
    done
done

./ms-stack.py $files
