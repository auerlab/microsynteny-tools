#!/bin/sh -e

make clean all
for gff in GFF/*.gff3; do
    printf "\n========================================================\n"
    printf "$gff\n"
    printf "========================================================\n"
    while read line; do
	printf "\n$line\n"
	gene=$(printf "$line\n" | cut -f 2);
	time ./msyn-hood --output-dir Hoods $gff $gene 200000
    done < DETF_refined.tsv
done
