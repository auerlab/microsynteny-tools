#!/bin/sh -e

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
    printf "Usage: $0 gene-name transcriptome|genome\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 2 ]; then
    usage
fi
gene=$1
db=$2

# Debug
# export PATH=~/Prog/Src/local/bin:$PATH

test -e $gene.fa || blt extract-seq GFF/Danio_rerio.GRCz11.105.chr.gff3 \
    Reference/Danio_rerio.GRCz11.dna.toplevel.fa.xz gene "Name=$gene;" \
    | grep -w -A 1 CDS > $gene.fa

# -task megablast (default) = hightly similar
# -task dc-megablast = discontinguous megablast (pretty similar_
# -task blastn = somewhat similar
# evalue < 0.01 considered good enough for homology
# https://www.metagenomics.wiki/tools/blast/evalue
# https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
case $db in
transcriptome)
    blastn -task blastn -db BLAST-DB/Homo-sapiens-transcriptome \
	-query $gene.fa -evalue 0.01 \
	-outfmt "6 qseqid pident evalue stitle" \
	> $gene-transcriptome.out 2> $gene-transcriptome.err
    printf "%-15s %10s %10s %-25s %s\n" \
	"Query" "Identity" "E-value" "Location" "Gene"
    awk -f blast-filter.awk $gene-transcriptome.out
    ;;

genome)
    if [ -e $gene-genome.out ]; then
	printf "$gene-genome.out already exists.  Rerun BLAST search? y/[n] "
	read run
    else
	run=y
    fi
    if [ 0$run = 0y ]; then
	# More than 2 threads backfires due to I/O
	time blastn -task dc-megablast -db BLAST-DB/Homo-sapiens-genome \
	    -query $gene.fa -outfmt "6 qseqid pident evalue stitle" \
	    -num_threads 2 \
	    > $gene-genome.out
    fi
    # blast-filter.awk will attempt to print gene name from $10 as
    # well, but it will just be blank
    printf "%-15s %10s %10s %-25s %s\n" \
	"Query" "Identity" "E-value" "Location"
    awk -f blast-filter.awk $gene-genome.out
    ;;

*)
    usage
    ;;

esac
