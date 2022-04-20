#!/bin/sh -e

while [ -z "$release" ]; do
    printf "Genome release? (e.g. 106) "
    read release
done

site="http://ftp.ensembl.org/pub/release-$release/gff3/danio_rerio"
curl_flags='--continue-at - --remote-name'
gff=Danio_rerio.GRCz11.$release.chr.gff3
if [ ! -e $gff ]; then
    echo $gff
    url=$site/$gff.gz
    curl $curl_flags $url
    gunzip $gff.gz
else
    printf "$(basename ${gff%.gz}) already exists.\n"
fi

awk '{ print $2 }' detf-refined.tsv | fgrep -v unnamed | sort --ignore-case \
    | uniq > ../Raw/de-tf-refined.txt

awk -F , 'length($2) > 0 && $1 ~ "ENSDART" { print $2 }' all-de.csv \
    | awk -f split-gene.awk | fgrep -v unnamed | sort --ignore-case \
    | uniq > ../Raw/all-de.txt

blt ensemblid2gene $gff detf-ids.txt \
    | awk '{ print $2 }' | fgrep -v unnamed | sort --ignore-case \
    | uniq > ../Raw/de-tf.txt

blt ensemblid2gene $gff expressed_NDE_transcripts.txt \
    | awk '{ print $2 }' | fgrep -v unnamed | sort --ignore-case \
    | uniq > ../Raw/all-non-de.txt

awk -F '[\t]' -f zf-nde-tfs.awk zf_NDE_TFs_Jake.csv | fgrep -v unnamed \
    | sort --ignore-case | uniq > non-de-tf-jake.txt
cut -d '|' -f 1 non-de-tf-jake.txt > ../Raw/non-de-tf.txt
