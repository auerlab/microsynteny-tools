#!/bin/sh -e

make clean all
while read line; do
    printf "\n$line\n"
    gene=$(printf "$line\n" | cut -f 2);
    time ./msyn-hood GFF/Danio_rerio.GRCz11.105.chr.gff3 $gene 200000
done < DETF_refined.tsv
