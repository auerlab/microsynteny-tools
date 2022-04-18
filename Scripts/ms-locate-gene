#!/bin/sh

##########################################################################
#   Description:
#       
#   History:
#   Date        Name        Modification
#   2022-03-01  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 gene-name GFF-file [GFF-file ...]\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# -lt 2 ]; then
    usage
fi
gene=$1
shift

for file in "$@"; do
    printf "\n$gene $file\n"
    awk -v gene=$gene '$3 == "gene" {
	split($9, a, ";");
	for (f in a) {
	    str = "Name=" gene;
	    #print a[f];
	    if ( tolower(a[f]) ~ tolower(str) ) {
		print $1, $4, a[f];
		break;
	    }
	}
    }' $file
done
