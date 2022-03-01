#!/bin/sh -e

for gene in $(cat Genes/all-de-names.txt); do
    ./common-to.sh $gene
done
