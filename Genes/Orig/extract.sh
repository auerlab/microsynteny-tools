#!/bin/sh -e

awk '{ print $2 }' detf-refined.tsv | fgrep -v unnamed | sort --ignore-case \
    | uniq > ../de-tf-refined.txt

awk -F , 'length($2) > 0 && $1 ~ "ENSDART" { print $2 }' all-de.csv \
    | awk -f split-gene.awk | fgrep -v unnamed | sort --ignore-case \
    | uniq > ../all-de.txt

blt ensemblid2gene \
    ../../GFF/Danio_rerio.GRCz11.105.chr.gff3 detf-ids.txt \
    | awk '{ print $2 }' | fgrep -v unnamed | sort --ignore-case \
    | uniq > ../de-tf.txt

blt ensemblid2gene \
    ../../GFF/Danio_rerio.GRCz11.105.chr.gff3 expressed_NDE_transcripts.txt \
    | awk '{ print $2 }' | fgrep -v unnamed | sort --ignore-case \
    | uniq > ../all-non-de.txt

awk -F '[\t]' -f zf-nde-tfs.awk zf_NDE_TFs_Jake.csv | fgrep -v unnamed \
    | sort --ignore-case | uniq > ../non-de-tf-jake.txt

