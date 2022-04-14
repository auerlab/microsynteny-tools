#!/bin/sh -e

if [ $0 != ./fetch-cdnas.sh ]; then
    printf "Must be run as ./fetch-cdnas.sh\n"
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

site="http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/cdna"

# No combined primary_assembly files for many species, so we would
# have to download all chromosomes individually and concatenate.
# Remove non-chromosomal sequences from toplevel to get primary_assembly.
for cdna in Homo_sapiens.GRCh38.cdna.all.fa
do
    base=$(basename $cdna)
    printf "===\n"
    if [ ! -e ${base%.gz} ]; then
	printf "Downloading ${base%.gz}\n"
	url=$site/$cdna
	$fetch $flags $url
    else
	printf "$(basename ${cdna%.gz}) already exists.\n"
    fi
done
