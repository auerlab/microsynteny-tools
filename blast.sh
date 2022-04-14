#!/bin/sh -ex

# Debug
# export PATH=~/Prog/Src/local/bin:$PATH

for gene in isl2b; do
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
done
