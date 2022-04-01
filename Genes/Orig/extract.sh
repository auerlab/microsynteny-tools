#!/bin/sh -e

awk '{ print $2 }' detf-refined.tsv | sort | uniq > ../de-tf-refined.txt

awk -F , 'length($2) > 0 && $1 ~ "ENSDART" { print $2 }' all-de.csv \
    | awk -f split-gene.awk | sort | uniq > ../all-de.txt

blt ensemblid2gene \
    ../../GFF/Danio_rerio.GRCz11.105.chr.gff3 detf-ids.txt \
    | awk '{ print $2 }' | sort | uniq > ../de-tf.txt

blt ensemblid2gene \
    ../../GFF/Danio_rerio.GRCz11.105.chr.gff3 expressed_NDE_transcripts.txt \
    | awk '{ print $2 }' | sort | uniq > ../all-non-de.txt

awk -f zf-nde-tfs.awk zf_NDE_TFs_Jake.csv | sort | uniq > ../non-de-tf-jake.txt

