#!/bin/sh -e

if [ $0 != ./fetch-cdnas.sh ]; then
    printf "Must be run as ./fetch-cdnas.sh\n"
    exit 1
fi

curl_flags='--continue-at - --remote-name'
site="http://ftp.ensembl.org/pub/release-$(../Utils/ensembl-release.sh)/fasta"

# No combined primary_assembly files for many species, so we would
# have to download all chromosomes individually and concatenate.
# Remove non-chromosomal sequences from toplevel to get primary_assembly.
for cdna in \
    homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa \
    mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa \
    rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa \
    danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa \
    takifugu_rubripes/cdna/Takifugu_rubripes.fTakRub1.2.cdna.all.fa \
    oryzias_latipes/cdna/Oryzias_latipes.ASM223467v1.cdna.all.fa
do
    base=$(basename $cdna)
    printf "===\n"
    printf "Downloading $base.gz\n"
    url=$site/$cdna.gz
    curl $curl_flags $url
done
