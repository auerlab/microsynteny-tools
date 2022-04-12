#!/bin/sh -e

##########################################################################
#   Description:
#       Extract genes and their neighbors from a GFF, saving the
#       neighborhood to another GFF.
#
#   History:
#   Date        Name        Modification
#   2022-02-04  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 [ms-extract-flags] gene-list.txt [file.gff3 ...]\n"
    ./ms-extract --help
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# -lt 1 ]; then
    usage
fi

if [ $# != 0 ]; then
    while [ 0$(echo $1 | cut -c 1) = 0'-' ]; do
	# All flags have one argument (for now)
	flags="$flags $1 $2"
	shift
	shift
    done
fi

gene_list=$1
shift

if [ $# = 0 ]; then
    files="GFF/Danio_rerio.GRCz11.105.chr.gff3"
    files="$files GFF/Oryzias_latipes.ASM223467v1.105.chr.gff3"
    files="$files GFF/Takifugu_rubripes.fTakRub1.2.105.chr.gff3"
    files="$files GFF/Mus_musculus.GRCm39.105.chr.gff3"
    files="$files GFF/Rattus_norvegicus.mRatBN7.2.105.chr.gff3"
    files="$files GFF/Homo_sapiens.GRCh38.105.chr.gff3"
else
    files="$@"
fi

printf "%s\n" "$flags"
printf "$files\n"

make clean install
mkdir -p Regions
for gff in $files; do
    printf "\n========================================================\n"
    printf "$gff\n"
    printf "========================================================\n"
    time ./ms-extract --output-dir Regions $flags $gff $gene_list
done
