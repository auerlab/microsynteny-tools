#!/bin/sh -e

if [ ! -e ../de-tf-refined.txt ]; then
    awk '{ print $2 }' detf-refined.tsv > ../de-tf-refined.txt
fi

if [ ! -e ../all-de.txt ]; then
    awk -F , 'length($2) > 0 && $1 ~ "ENSDART" { print $2 }' all-de.csv \
	| awk -f split-gene.awk > ../all-de.txt
fi

if [ ! -e ../de-tf.txt ]; then
    blt ensemblid2gene \
	../../GFF/Danio_rerio.GRCz11.105.chr.gff3 detf-ids.txt \
	| awk '{ print $2 }' | sort | uniq > ../de-tf.txt
fi

if [ ! -e ../all-non-de.txt ]; then
    blt ensemblid2gene \
	../../GFF/Danio_rerio.GRCz11.105.chr.gff3 expressed_NDE_transcripts.txt \
	| awk '{ print $2 }' | sort | uniq > ../all-non-de.txt
fi

if [ ! -e ../non-de-tf.txt ]; then
    awk '{ print $1 }' zf_NDE_TFs_Jake.csv > ../non-de-tf.txt
fi
