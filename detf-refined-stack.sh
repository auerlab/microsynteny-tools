#!/bin/sh -e

for gene in $(cat Genes/detf-refined-names.txt); do
    ./stack.sh $gene
done
