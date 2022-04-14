#!/bin/sh -e

if [ $0 != ./fetch-refs.sh ]; then
    printf "Must be run as ./fetch-refs.sh\n"
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

site="http://ftp.ensembl.org/pub/release-105/fasta"
# No combined primary_assembly files for many species, so we would
# have to download all chromosomes individually and concatenate.
# Remove non-chromosomal sequences from toplevel to get primary_assembly.
ref_type="toplevel"
for ref in \
    danio_rerio/dna/Danio_rerio.GRCz11.dna.$ref_type.fa.gz \
    oryzias_latipes/dna/Oryzias_latipes.ASM223467v1.dna.$ref_type.fa.gz \
    takifugu_rubripes/dna/Takifugu_rubripes.fTakRub1.2.dna.$ref_type.fa.gz \
    mus_musculus/dna/Mus_musculus.GRCm39.dna.$ref_type.fa.gz \
    rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.$ref_type.fa.gz \
    homo_sapiens/dna/Homo_sapiens.GRCh38.dna.$ref_type.fa.gz
do
    base=$(basename $ref)
    printf "===\n"
    if [ ! -e ${base%.gz} ]; then
	printf "Downloading ${base%.gz}\n"
	url=$site/$ref
	$fetch $flags $url
    else
	printf "$(basename ${ref%.gz}) already exists.\n"
    fi
done
