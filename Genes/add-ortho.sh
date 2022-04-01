#!/bin/sh -e

cd ..
make clean all
cd Genes
cut -d '|' -f 1 non-de-tf-jake.txt > non-de-tf.txt
for file in all-de.txt de-tf-refined.txt non-de-tf.txt \
	    all-non-de.txt de-tf.txt; do
    echo $file
    ../add-ortho $file biomart-orthologs-sorted.tsv.xz \
	> ${file%.txt}-ortho.txt
done
