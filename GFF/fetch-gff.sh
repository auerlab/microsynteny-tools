#!/bin/sh -e

if [ $0 != ./fetch-gff.sh ]; then
    printf "Must be run as ./fetch-gff.sh\n"
    exit 1
fi

for fetch in curl wget fetch; do
    which $fetch > /dev/null 2>&1 && break
done
case $fetch in
curl)
    flags='--continue-at - --remote-name'
    ;;
wget)
    flags='--continue'
    ;;
*)
    ;;
esac

site='http://ftp.ensembl.org/pub/release-105/gff3'
for gff in \
    takifugu_rubripes/Takifugu_rubripes.fTakRub1.2.105.chr.gff3.gz \
    danio_rerio/Danio_rerio.GRCz11.105.chr.gff3.gz \
    xenopus_tropicalis/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.105.chr.gff3.gz \
    mus_musculus/Mus_musculus.GRCm39.105.chr.gff3.gz; do
    if [ ! -e $(basename ${gff%.gz}) ]; then
	url=$site/$gff
	$fetch $flags $url
	gunzip $(basename $gff)
    else
	printf "$(basename ${gff%.gz}) already exists.\n"
    fi
done
