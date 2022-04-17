#!/bin/sh -e

if [ $0 != ./fetch-refs.sh ]; then
    printf "Must be run as ./fetch-refs.sh\n"
    exit 1
fi

curl_flags='--continue-at - --remote-name'
site="http://ftp.ensembl.org/pub/release-$(../Utils/ensembl-release.sh)/fasta"

# No combined primary_assembly files for many species, so we would
# have to download all chromosomes individually and concatenate.
# Remove non-chromosomal sequences from toplevel to get primary_assembly.
ref_type="toplevel"
for ref in \
    danio_rerio/dna/Danio_rerio.GRCz11.dna.$ref_type.fa \
    oryzias_latipes/dna/Oryzias_latipes.ASM223467v1.dna.$ref_type.fa \
    takifugu_rubripes/dna/Takifugu_rubripes.fTakRub1.2.dna.$ref_type.fa \
    mus_musculus/dna/Mus_musculus.GRCm39.dna.$ref_type.fa \
    rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.$ref_type.fa \
    homo_sapiens/dna/Homo_sapiens.GRCh38.dna.$ref_type.fa
do
    base=$(basename $ref)
    printf "===\n"
    printf "Downloading $base...\n"
    url=$site/$ref.gz
    curl $curl_flags $url
done
