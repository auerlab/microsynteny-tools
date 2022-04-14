#!/bin/sh -ex

# Debug
export PATH=~/Prog/Src/local/bin:$PATH

blt extract-seq GFF/Danio_rerio.GRCz11.105.chr.gff3 \
    Reference/Danio_rerio.GRCz11.dna.toplevel.fa.xz gene 'Name=isl2b;' \
    | grep -w -A 1 CDS > query.fa
blastn -db BLAST-DB/Homo-sapiens -query query.fa
