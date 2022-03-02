#!/bin/sh -e

for gene in $(cat Genes/detf-refined-names.txt); do
    ./intersect.sh $gene
done
