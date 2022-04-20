#!/bin/sh -e

if [ ! -e biomart-orthologs.tsv.xz ]; then
    cat << EOM

Missing biomart-orthologs.tsv.xz. Go to:

https://www.ensembl.org/biomart/martview/

Choose database:    Ensembl Genes
Choose dataset:     Zebrafish genes
Under "Gene":
    Clear all attributes
Under "Homologs:"
    Select Homologs/GENE/Gene name
    For each species (Zebrafish, Medaka, Fugu, Mouse, Rat, Human),
    select Homologs/Species/Gene name
    Site says up to 6 orthologs, but consistently fails for 5 if trying
    to immediately save to a file.
    Choose "Compressed web file (notify by email)" and expect to wait a
    few hours.

    Copy downloaded file to biomart-orthologs.tsv and compress with "xz".

EOM
fi

if [ ! -e biomart-orthologs-sorted.tsv.xz ]; then
    printf "Sorting...\n"
    ./sort-biomart.sh
fi

cd ..
./cave-man-install.sh
cd Genes
for file in Raw/*.txt; do
    echo $file
    ms-add-ortho $file biomart-orthologs-sorted.tsv.xz \
	> $(basename ${file%.txt}-ortho.txt)
done
