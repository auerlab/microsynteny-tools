#!/bin/sh -e

##########################################################################
#   Description:
#       
#   History:
#   Date        Name        Modification
#   2022-03     Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 [ms-stack.py flags] adjacent-genes max-nt gene-list.txt\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# -lt 3 ]; then
    usage
fi
while [ $(echo $1 | cut -c 1,1) = '-' ]; do
    flags="$flags $1"
    shift
done

adjacent_genes=$1
max_nt=$2
gene_file="$3"

for gene in $(cat $gene_file); do
    ./stack.sh $flags --show-gene-lens $adjacent_genes $max_nt \
	$(echo $gene | tr '|' ' ')
done
