#!/bin/sh -e

##########################################################################
#   Description:
#       This script is not currently useful.  See readme and use
#       biomart-generated ortholog lists.
#       
#   History:
#   Date        Name        Modification
#   2022-04-01  Jason Bacon Begin
##########################################################################

# https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/
# README_ensembl:   TaxID -> species name
# gene_info:        TaxID, GeneID, synonyms, etc.
# gene_orthologs:   TaxID, Gene_ID -> TaxID, GeneID

for file in gene_info.gz gene_orthologs.gz README_ensembl; do
    if [ ! -e $file ]; then
	curl -O https://ftp.ncbi.nlm.nih.gov/gene/DATA/$file
    else
	printf "$file already exists.\n"
    fi
done

# README_ensembl
# TaxID     Species
# 7955      Danio Rerio (Zebrafish)
# 8090      Oryzias latipes (Medaka)
# 31033     Takifugu rubripes (Puffer)
# 10090     Mus Musculus
# 10116     Rattus Norvegicus
# 9606      Home sapiens

#zmore gene_orthologs.gz
#zmore gene_info.gz

goi=isl2b
#goi=sox11a

printf "Getting TaxId and GeneID...\n"
ids=$(zcat gene_info.gz | awk -v goi=$goi '$1 == 7955 && $3 == goi { print $1, $2 }')
taxid=$(echo $ids | cut -d ' ' -f 1)
geneid=$(echo $ids | cut -d ' ' -f 2)
echo $taxid $geneid

printf "Listing orthologs...\n"
zcat gene_orthologs.gz \
    | awk -v taxid=$taxid -v geneid=$geneid '$1 == taxid && $2 == geneid' \
    > $goi-ortho.txt
cat $goi-ortho.txt

rm -f interesting-ortho.txt
for taxid in 8090 31033 10090 10116 9606; do
    awk -v taxid=$taxid '$4 == taxid' $goi-ortho.txt >> interesting-ortho.txt
done

printf "Ortholog gene names...\n"
while read line; do
    taxid=$(echo $line | cut -d ' ' -f 4)
    geneid=$(echo $line | cut -d ' ' -f 5)
    echo $taxid $geneid
    zcat gene_info.gz \
	| awk -v taxid=$taxid -v geneid=$geneid '$1 == taxid && $2 == geneid'
done < interesting-ortho.txt
