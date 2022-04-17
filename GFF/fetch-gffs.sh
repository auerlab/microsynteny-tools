#!/bin/sh -e

if [ $0 != ./fetch-gffs.sh ]; then
    printf "Must be run as ./fetch-gffs.sh\n"
    exit 1
fi

#    gadus_morhua/Gadus_morhua.gadMor3.0.105.chr.gff3.gz \
#    salmo_salar/Salmo_salar.ICSASG_v2.105.chr.gff3.gz \
#    gouania_willdenowi/Gouania_willdenowi.fGouWil2.1.105.chr.gff3.gz \
#    salmo_trutta/Salmo_trutta.fSalTru1.1.105.chr.gff3.gz \
#    cottoperca_gobio/Cottoperca_gobio.fCotGob3.1.105.chr.gff3.gz \
#    ictalurus_punctatus/Ictalurus_punctatus.IpCoco_1.2.105.chr.gff3.gz
#    xenopus_tropicalis/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.105.chr.gff3 \
#    gallus_gallus/Gallus_gallus.GRCg6a.105.chr.gff3 \

site='http://ftp.ensembl.org/pub/release-105/gff3'
curl_flags='--continue-at - --remote-name'
for gff in \
    danio_rerio/Danio_rerio.GRCz11.105.chr.gff3 \
    oryzias_latipes/Oryzias_latipes.ASM223467v1.105.chr.gff3 \
    takifugu_rubripes/Takifugu_rubripes.fTakRub1.2.105.chr.gff3 \
    mus_musculus/Mus_musculus.GRCm39.105.chr.gff3 \
    rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.105.chr.gff3 \
    homo_sapiens/Homo_sapiens.GRCh38.105.chr.gff3
do
    base=$(basename $gff)
    if [ ! -e $base ]; then
	echo $gff
	url=$site/$gff.gz
	curl $curl_flags $url
	gunzip $base.gz
    else
	printf "$(basename ${gff%.gz}) already exists.\n"
    fi
done
