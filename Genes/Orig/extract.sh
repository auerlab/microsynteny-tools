#!/bin/sh -e

awk '{ print $2 }' detf-refined.tsv > ../detf-refined.txt

awk -F , 'length($2) > 0 && $1 ~ "ENSDART" { print $2 }' all-de.csv \
    | awk -f split-gene.awk > ../all-de.txt

blt ensemblid2gene \
    ../../GFF/Danio_rerio.GRCz11.105.chr.gff3 detf-ids.txt \
    | awk '{ print $2 }' | sort | uniq > ../detf.txt

blt ensemblid2gene \
    ../../GFF/Danio_rerio.GRCz11.105.chr.gff3 expressed_NDE_transcripts.txt \
    | awk '{ print $2 }' | sort | uniq > ../non-de.txt
