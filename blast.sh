#!/bin/sh -ex

##########################################################################
#   Description:
#       BLAST a Danio gene against other transcriptomes/genomes
#       
#   History:
#   Date        Name        Modification
#   2022-04-15  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 gene-name\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
gene=$1

# Debug
# export PATH=~/Prog/Src/local/bin:$PATH

test -e $gene.fa || blt extract-seq GFF/Danio_rerio.GRCz11.105.chr.gff3 \
    Reference/Danio_rerio.GRCz11.dna.toplevel.fa.xz gene "Name=$gene;" \
    | grep -w -A 1 CDS > $gene.fa

# -task megablast (default) = hightly similar
# -task dc-megablast = discontinguous megablast (pretty similar_
# -task blastn = somewhat similar
blastn -task dc-megablast -db BLAST-DB/Homo-sapiens-transcriptome \
    -query $gene.fa > $gene-transcriptome.out
more $gene-transcriptome.out

blastn -task dc-megablast -db BLAST-DB/Homo-sapiens-genome \
    -query $gene.fa > $gene-genome.out
more $gene-genome.out
