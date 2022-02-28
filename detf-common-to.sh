#!/bin/sh -e

for gene in $(cat Genes/detf-names.txt); do
    ./common-to.sh $gene
done
